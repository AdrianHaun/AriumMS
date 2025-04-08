classdef RawData
    %% RawData superclass object for AriumMS file processing
    % Processing includes:
    % - Loading of MS data files
    % - Pre processing
    % - ROI search
    % - Post processing
    % - Feature generation

    properties
        %% Processing Parameters
        fileName    (:,1) string
        dataFile    (:,1) string
        blankFile   (:,1) string
        groupName   (1,1) string
        % Main Processing Options
        useBlankSubtraction   (1,1) logical = false
        useSmoothing          (1,1) logical = false
        useBaselineCorrection (1,1) logical = false
        useIsotopeFilter      (1,1) logical = false
        useAdductFilter       (1,1) logical = false
        useContaminantFilter  (1,1) logical = false
        useInternalStandard   (1,1) logical = false
        useMassAlign          (1,1) logical = false
        usePeakAlign          (1,1) logical = false
        useScaling            (1,1) logical = false
        % ROI parameter
        roiThreshold            (1,1) double {mustBeInteger,mustBePositive} = 5000
        withinFileMassTolerance (1,1) double {mustBePositive} = 0.01
        withinFileMassUnit      (1,1) string {mustBeMember(withinFileMassUnit,["Da","ppm"])} = "Da"
        roiMinOccurrence         (1,1) double {mustBeInteger,mustBePositive} = 20
        measurementStartTime    (1,1) double {mustBeFinite} = 0
        measurementEndTime      (1,1) double {mustBeFinite} = 1
        % Baseline Correction Parameters
        baselineWindowSize          (1,1) double {mustBeFinite,mustBePositive} = 20
        baselineStepSize            (1,1) double {mustBeFinite,mustBePositive} = 20
        baselineRegressionMethod    (1,1) string {mustBeMember(baselineRegressionMethod,["pchip","linear","spline"])} = "pchip"
        baselineEstimationMethod    (1,1) string {mustBeMember(baselineEstimationMethod,["quantile","em"])} = "em"
        baselineSmoothMethod        (1,1) string {mustBeMember(baselineSmoothMethod,["none","lowess","loess"])} = "none"
        baselineQuantil             (1,1) double {mustBeInRange(baselineQuantil,0,1)} = 0.1
        %Golay Parameters
        smoothingFrameSize  (1,1) double {mustBeInteger,mustBePositive} = 20
        smoothingDegree     (1,1) double {mustBeInteger,mustBePositive} = 2
        % Internal Standard Data
        nInternalStandard     (1,1) double {mustBeInteger,mustBePositive} = 1
        internalStandardData  (:,3) double
        useISMassCorrection   (1,1) logical = false
        applyISto             (1,1) string {mustBeMember(applyISto,["S&B","SOnly"])} = "SOnly"
        internalStandardOrder (1,1) string {mustBeMember(internalStandardOrder,["BlankIS","ISBlank"])} = "ISBlank"
        %Adduct Parameters
        scanPolarity        (1,1) string {mustBeMember(scanPolarity,["positive","negative","both"])} = "positive"
        minCosineSimilarity (1,1) double {mustBeInRange(minCosineSimilarity,0,1)} = 0.85
        adductSelectedPos   (30,1) logical = false        %Structure: 1:12 Single Charged, 13:18 Dimers, 19:26 DoubleCharged, 27:30 TripleCharged
        adductSelectedNeg   (16,1) logical = false       %Structure: 1:10 SingleCharged, 11:14 Dimers, 15 DoubleCharged, 16 TripleCharged
        neutralSelectedSmol (18,1) logical = false
        neutralSelectedCon  (17,1) logical = false
        % Mass Alignment Parameters
        massAlignmentEstimMethod        (1,1) string {mustBeMember(massAlignmentEstimMethod,["histogram","regression"])} = "regression"
        massAlignmentCorrectionMethod   (1,1) string {mustBeMember(massAlignmentCorrectionMethod,["nearest-neighbor","shortest-path"])} = "nearest-neighbor"
        massAlignmentQuantil            (1,1) double {mustBeInRange(massAlignmentQuantil,0,1)} = 0.99
        % Peak Align Parameters
        peakAlignmentMaxShiftNegative   (1,1) double {mustBeFinite,mustBePositive} = 20
        peakAlignmentMaxShiftPositive   (1,1) double {mustBeFinite,mustBePositive} = 20
        peakAlignmentPulseWidth         (1,1) double {mustBeFinite,mustBePositive} = 2
        peakAlignmentWindowSizeRatio    (1,1) double {mustBePositive} = 2.5
        peakAlignmentSearchSpace        (1,1) string {mustBeMember(peakAlignmentSearchSpace,["regular","latin"])} = "regular"
        peakAlignmentIteration         (1,1) double {mustBeInteger,mustBePositive} = 5
        peakAlignmentGridSteps          (1,1) double {mustBeInteger,mustBePositive} = 20
        % Scaling
        groupScale          (1,1) double {mustBePositive} = 1
        sampleScale         (:,1) double {mustBePositive} = 1
        % Integration and Filter
        evaluationParameter         (1,1) string {mustBeMember(evaluationParameter,["Area","Height"])} = "Area"
        peakMinWidth                (1,1) double {mustBePositive} = 2
        peakMaxWidth                (1,1) double {mustBePositive} = 45
        minSignal2Noise             (1,1) double {mustBePositive} = 3
        minOccurence                (1,1) double {mustBeInRange(minOccurence,0,1)} = 0.5
        betweenFileMassTolerance    (1,1) double {mustBeFinite} = 0.001
        betweenFileMassUnit         (1,1) string {mustBeMember(betweenFileMassUnit,["Da","ppm"])} = "Da"
        peakTimeTolerance           (1,1) double {mustBeFinite} = 5
        useEntropyFilter            (1,1) logical = false
        entropyFilterStrength       (1,1) string {mustBeMember(entropyFilterStrength,["strict","medium","lax"])} = "medium"
        %% DataStorage
        RawDataFile         string
        RawDataFileObj      (1,1)
        TempDataFile        string
        TempDataFileObj     (1,1)
        ROIDataFile         string
        ROIDataFileObj      (1,1)
        % processing variables
        nScan              (:,1) double {mustBeInteger,mustBePositive}
        nScanPadded        (:,1) double {mustBeInteger,mustBePositive}
        %FileInfos
        DataInfo            (1,:) struct
        scanFrequencySecond  (1,1) double
        %IS Data
        interalStandardIntensity    (:,:) double
        internalStandardTime        (:,:) double
        internalStandardMass        (1,:) double
        internalStandardMassFound   (1,:) double
        internalStandardMassDelta   (:,:) double
        mzCorrectionFcn             (1,1)
        %% Plot
        mainWindow          matlab.ui.Figure

        %testing variables
        Output

    end

    methods
        function obj = RawData(groupName,appWindow)
            % Construct an instance of this class and generate semi
            % permanant storage files
            if isgraphics(appWindow)
                obj.mainWindow = appWindow;
            end
            obj.groupName = groupName;
            obj.RawDataFile = tempname +".mat";
            obj.ROIDataFile = tempname +".mat";
            obj = obj.initializeStorageFile;            
        end
        %% Handling MS files
        function obj = dataFileCheck(obj)
            % Check Data, number of Scans, Start/End Times
            %Check minimum number of scans
            fileArray = [obj.dataFile;obj.blankFile];
            %remove empty
            idx = cellfun(@isempty,fileArray);
            fileArray(idx) = [];
            %preallocation
            retentionTimes = cell(length(fileArray),1);
            totalIonChromatogram = cell(length(fileArray),1);
            basePeakChromatogram = cell(length(fileArray),1);
            polarityCells = cell(length(fileArray),1);
            FileInfo =  struct('numberOfScansMS1',[],...
                'numberOfScansMSn',[],...
                'startTime',[],...
                'endTime',[],...
                'scanFrequenceMS1',[],...
                'scanFrequenceMS2',[]);

            parfor iFile = 1:numel(fileArray)
                %check file-type
                test = strsplit(fileArray(iFile),'.');
                test = test(end);
                switch test
                    case "mzML"
                        [FileInfo(iFile),retentionTimes{iFile},totalIonChromatogram{iFile},basePeakChromatogram{iFile},polarityCells{iFile}] = mzMLinfo(fileArray{iFile});
                    case "mzXML"
                        [FileInfo(iFile),retentionTimes{iFile},totalIonChromatogram{iFile},basePeakChromatogram{iFile},polarityCells{iFile}] = mzXMLinfo(fileArray{iFile});
                    case "CDF"
                        [FileInfo(iFile),retentionTimes{iFile},totalIonChromatogram{iFile},basePeakChromatogram{iFile},polarityCells{iFile}] = mzCDFinfo(fileArray{iFile});
                end
            end
            obj.RawDataFileObj.previewTICs = totalIonChromatogram;
            obj.RawDataFileObj.previewBPCs = basePeakChromatogram;
            obj.RawDataFileObj.previewTimes = retentionTimes;

            polarity = vertcat(polarityCells{:});

            test = strcmp(polarity,"+");
            if all(test)
                polarity = "positive";
            elseif all(~test)
                polarity = "negative";
            else
                polarity = "both";
            end
            obj.scanPolarity = polarity;

            % calculate Scan Frequency [Hz]
            scanFrequencyHertz = [FileInfo.numberOfScansMS1]./([FileInfo.endTime]-[FileInfo.startTime]);
            scanFrequencyHertz = num2cell(scanFrequencyHertz);
            [FileInfo.scanFrequencyMS1] = scanFrequencyHertz{:};
            scanFrequencyHertz = [FileInfo.numberOfScansMSn]./([FileInfo.endTime]-[FileInfo.startTime]);
            scanFrequencyHertz = num2cell(scanFrequencyHertz);
            [FileInfo.scanFrequencyMSn] = scanFrequencyHertz{:};
            %store data
            obj.DataInfo = FileInfo;
            obj.measurementStartTime = round(min([FileInfo.startTime]),1);
            obj.measurementEndTime = round(max([FileInfo.endTime]),1);
        end

        function obj = readData(obj,dataFile,separationType)
            nFile = size(dataFile,1);
            %preallocation
            DataMS1 = cell(nFile,1);
            DataMS2 = cell(nFile,1);
            parfor iFile = 1:nFile
                %filetype check
                fileType = strsplit(dataFile(iFile),'.');
                fileType = fileType(end);
                switch fileType
                    case "mzML"
                        [DataMS1{iFile,1},DataMS2{iFile,1}] = readmzML_MSandMS2(dataFile{iFile});
                    case "mzXML"
                        [DataMS1{iFile,1},DataMS2{iFile,1}] = readmzXML_MSandMS2(dataFile{iFile});
                    case "CDF"
                        [DataMS1{iFile,1},DataMS2{iFile,1}] =  readmzCDF(dataFile{iFile});
                end
            end
            
            % pre-process raw data
            parfor iFile = 1:nFile
                switch separationType
                    case "GC"
                        %centroid profile data
                        DataMS1{iFile,1}.centroidDataMS1 = centroidScans(DataMS1{iFile,1}.profileDataMS1);
                        % clean scans
                        DataMS1{iFile,1}.centroidDataMS1 = cleanScans(DataMS1{iFile,1}.centroidDataMS1);
                        DataMS1{iFile,1}.profileDataMS1 = cleanScans(DataMS1{iFile,1}.profileDataMS1);

                    otherwise
                        %convert MS1 and MS2 precursor to molecular mass
                        DataMS1{iFile,1}.profileDataMS1 = convertScans2MolecularMass(DataMS1{iFile,1}.profileDataMS1,DataMS1{iFile,1}.polarityMS1);

                        precursors = DataMS2{iFile,1}.precursorMass;
                        modifier = ones(size(precursors))*1.007825;
                        idx = DataMS2{iFile,1}.polarityMS2 == "+";
                        modifier(idx) = modifier(idx)*-1;
                        DataMS2{iFile,1}.precursorMassCorrected = precursors + modifier;

                        %centroid profile data
                        DataMS1{iFile,1}.centroidDataMS1 = centroidScans(DataMS1{iFile,1}.profileDataMS1);
                        %clean scans
                        DataMS1{iFile,1}.centroidDataMS1 = cleanScans(DataMS1{iFile,1}.centroidDataMS1);
                        DataMS1{iFile,1}.profileDataMS1 = cleanScans(DataMS1{iFile,1}.profileDataMS1);

                        %compress MS2 data and store
                        DataMS2{iFile,1}.centroidDataMS2 = centroidScans(DataMS2{iFile,1}.profileDataMS2);
                        DataMS2{iFile,1}.centroidDataMS2 = normalizeScans(DataMS2{iFile,1}.centroidDataMS2);
                        DataMS2{iFile,1}.centroidDataMS2 = cleanScans(DataMS2{iFile,1}.centroidDataMS2);
                end
            end
            %store data
            %MS1 data
            DataMS1 = vertcat(DataMS1{:});
            obj.RawDataFileObj.profileDataMS1 = {DataMS1.profileDataMS1}';
            obj.RawDataFileObj.centroidedDataMS1 = {DataMS1.centroidDataMS1}';
            obj.RawDataFileObj.timeDataMS1 = {DataMS1.timeDataMS1}';
            obj.RawDataFileObj.polarityMS1 = {DataMS1.polarityMS1}';

            %check for empty MS2 data
            DataMS2 = vertcat(DataMS2{:});
            if ~isscalar(vertcat(DataMS2(:).profileDataMS2))
                obj.RawDataFileObj.profileDataMS2 = {DataMS2.profileDataMS2}';
                obj.RawDataFileObj.centroidedDataMS2 = {DataMS2.centroidDataMS2}';
                obj.RawDataFileObj.timeDataMS2 = {DataMS2.timeDataMS2}';
                obj.RawDataFileObj.polarityMS2 = {DataMS2.polarityMS2}';
                obj.RawDataFileObj.precursorMass = {DataMS2.precursorMass}';
                obj.RawDataFileObj.molecularPrecursorMass = {DataMS2.precursorMassCorrected}';
                obj.RawDataFileObj.fragmentationEnergy = {DataMS2.fragmentationEnergy}';
                obj.RawDataFileObj.fragmentationType = {DataMS2.fragmentationType}';
            end
        end

        %% Data Processing
        function obj = averageBlankFiles(obj)
            %% averageBlankFiles calculates an average file from all blank files
            % Takes all blank files in the TempDataFile, calculates the
            % average and replaces blank files with the average blank.

            %gather data
            N_BLANK = numel(obj.blankFile);
            roiCell = obj.TempDataFileObj.roiCells;
            timeCell = obj.TempDataFileObj.timeCells;
            maxScan = max(obj.nScan);
            blankFiles = vertcat(roiCell{end-N_BLANK+1:end});
            blankFiles = reshape(blankFiles,maxScan,size(blankFiles,2),N_BLANK);
            blankFiles = mean(blankFiles,3);
            blankTimes = horzcat(timeCell{end-N_BLANK+1:end});
            %replace 0 with NaN then ignore NaN in median calculation
            blankTimes(blankTimes == 0) = NaN;
            blankTimes = median(blankTimes,2,"omitnan");
            blankTimes(isnan(blankTimes)) = 0;
            %replace blank data with average blank
            roiCell(end-N_BLANK+1:end) = [];
            timeCell(end-N_BLANK+1:end) = [];
            roiCell{end+1} = blankFiles;
            timeCell{end+1} = blankTimes;
            %store on disk
            obj.TempDataFileObj.roiCells = roiCell;
            obj.TempDataFileObj.timeCells = timeCell;
        end

        function obj = removeContaminants(obj)
            %% removeContaminants identifies contaminant masses and removes them
            % Identifies contaminants based on the UWPR common ESI
            % contaminants list. Identified mass is removed from found ROI
            % list.

            POLARITY = obj.scanPolarity;
            switch POLARITY
                case "positive"
                    contaminants = load("MassListData.mat","ContaminantsPos");
                    contaminants = contaminants.ContaminantsPos;
                case "negative"
                    contaminants = load("MassListData.mat","ContaminantsNeg");
                    contaminants = contaminants.ContaminantsNeg;
                otherwise
                    contaminants = load("MassListData.mat","ContaminantsPos","ContaminantsNeg");
                    contaminants = unique([contaminants.ContaminantsPos;contaminants.ContaminantsNeg]);
            end
            %calculate possible Contaminants
            switch obj.withinFileMassUnit
                case "Da"
                    isContaminant = abs(contaminants-obj.TempDataFileObj.ROImzVec) <= obj.withinFileMassTolerance;
                case "ppm"
                    isContaminant = abs(contaminants-obj.TempDataFileObj.ROImzVec)./obj.TempDataFileObj.ROImzVec*10^6 <= obj.withinFileMassTolerance;
            end
            isContaminant = any(isContaminant,1);
            % remove contaminant columns from ROI mz list and MSroi matrices
            obj.TempDataFileObj.ROImzVec(isContaminant) = [];
            tempCell = obj.TempDataFileObj.ROICells;
            for iTempCell = 1:size(obj.TempDataFileObj.ROICells,1)
                tempCell{iTempCell,1}(:,isContaminant)=[];
            end
            obj.TempDataFileObj.ROICells = tempCell;
        end

        function obj = internalStandardNormalization(obj)
            %% MUST BE UPDATED TO NEW PROCESSING
            %Gather relevant matrices
            if obj.applyISto == "S&B"
                Data = mat2cell(obj.TempDataFileObj.ROIMat,obj.nScanPadded);
                Data{end+1} = obj.TempDataFileObj.ROIMatBLK;
            else
                Data = mat2cell(obj.TempDataFileObj.ROIMat,obj.nScanPadded);
            end

            timeVectors = mat2cell(obj.TempDataFileObj.timeVec,obj.nScanPadded);

            % normalize Intensities
            if isscalar(obj.interalStandardIntensity) %only one IS
                value = obj.interalStandardIntensity;
                parfor id = 1:size(Data,1)
                    Data{id} = Data{id}./value;
                end

            else % multiple IS - Retention time dependent IS normalization
                retentionTimes = obj.internalStandardTime;
                intensities = obj.interalStandardIntensity;
                % Set up fittype and options.
                ft = 'pchipinterp';
                opts = fitoptions( 'Method', 'PchipInterpolant' );
                opts.ExtrapolationMethod = 'nearest';
                parfor id = 1:size(Data,1)
                    %% fit correction function
                    [xData, yData] = prepareCurveData(retentionTimes, intensities(id,:));
                    % Fit model to data.
                    intFcn = fit(xData,yData,ft,opts);
                    correctionVector = intFcn(timeVectors{id});
                    Data{id} = Data{id}./correctionVector;
                end
            end

            %store corrected Matrices back into object
            if obj.applyISto == "S&B"
                obj.TempDataFileObj.ROIMatBLK = Data{end};
                Data(end) = [];
            end
            obj.TempDataFileObj.ROIMat = vertcat(Data{:});
        end

        function obj = massCorrectionByInternalStandard(obj)
            %% massCorrectionByInternalStandard corrects found ROI masses
            % Takes the difference in found Internal Standard mass and
            % specified mass and shifts ROI masses accordingly.
            % Uses a fixed value for a single IS
            % Uses Pchip interpolation and nearest point extrapolation for
            % multiple IS

            if isscalar(obj.internalStandardMassDelta) % single IS constant correction
                obj.TempDataFileObj.ROImzVec = obj.TempDataFileObj.ROImzVec - obj.internalStandardMassDelta;

            else  % multiple IS - m/z dependent correction
                %% fit correction function
                [massData, massDelta] = prepareCurveData(obj.internalStandardMass, obj.internalStandardMassDelta);
                % Set up fittype and options.
                FUNCTION_TYPE = 'pchipinterp';
                Options = fitoptions( 'Method', 'PchipInterpolant' );
                Options.ExtrapolationMethod = 'nearest';
                % Fit model to data.
                obj.mzCorrectionFcn = fit(massData,massDelta,FUNCTION_TYPE,Options);

                % build mass correction vector and subtract from ROI masses
                mzCorrectionVector = obj.mzCorrectionFcn(obj.TempDataFileObj.ROImzVec);
                obj.TempDataFileObj.ROImzVec = obj.TempDataFileObj.ROImzVec - mzCorrectionVector;
            end
        end

        function [valuesFiltered,obj] = filterAdducts(obj,IntegrationResults)
            %% MUST BE UPDATED TO NEW PROCESSING
            % AdductFilterAlgo Filters Adduct Peaks from Internal AriumMS integration results
            %   Calculates possible non Adduct (Base) m/z for each
            %   extracted m/z, then finds matching masses in original list.
            %   Peak lists of Adduct and Base m/z are then compared to have
            %   peaks with the same RT (within 2 sec). Matching peaks
            %   shapes are compared using Cosine Similarity (>=0.85
            %   default)

            mzTolVal = obj.betweenFileMassTolerance;

            polarity = app.Data(CurrentTab).RawDataFileObj.polarity;
            polarity = vertcat(polarity{:});
            test = strcmp(polarity,"+");
            if all(test)
                polarity = "positive";
            elseif all(~test)
                polarity = "negative";
            else
                polarity = "both";
            end

            switch polarity
                case "positive"
                    Rules = load("MassListData.mat","AddPosRules","NLossRules");
                    Rules = [Rules.AddPosRules(obj.adductSelectedPos);Rules.NLossRules([obj.neutralSelectedSmol;obj.neutralSelectedCon])];
                case "negative"
                    Rules = load("MassListData.mat","AddNegRules","NLossRules");
                    Rules = [Rules.AddNegRules(obj.adductSelectedNeg);Rules.NLossRules([obj.neutralSelectedSmol;obj.neutralSelectedCon])];
                case "both"
                    Rules = load("MassListData.mat","AddNegRules","AddNegRules","NLossRules");
                    Rules = [Rules.AddPosRules(obj.adductSelectedPos);Rules.AddNegRules(obj.adductSelectedNeg);Rules.NLossRules([obj.neutralSelectedSmol;obj.neutralSelectedCon])];
            end
            %Preparation
            minCosSim = obj.minCosineSimilarity;
            mzVec = obj.TempDataFileObj.ROImzVec;
            EICMat = obj.TempDataFileObj.ROIMat;
            RTs=cellfun(@(X) X(:,2),IntegrationResults(3,:),'UniformOutput',false); %extract Retentiontimes
            ranges=cellfun(@(X) X(:,2:3),IntegrationResults(2,:),'UniformOutput',false); %extract Peak ranges

            %build possible Adduct mass list
            PossibleBaseMZ=zeros(length(Rules),length(mzVec));
            parfor n=1:length(Rules)
                PossibleBaseMZ(n,:)=Rules{n}(mzVec);
            end

            %compare case masses to initial mz list
            switch obj.betweenFileMassUnit
                case "ppm"
                    PossibleBaseMZindex = cell(size(PossibleBaseMZ));
                    for n=1:size(mzVec,2)
                        [~,PossibleBaseMZindex(:,n)]=ismembertol(PossibleBaseMZ(:,n),mzVec,mzTolVal,'DataScale',mzVec(n)/10^6,'OutputAllIndices',true);
                    end
                case "Da"
                    [~,PossibleBaseMZindex]=ismembertol(PossibleBaseMZ,mzVec,mzTolVal,'DataScale',1,'OutputAllIndices',true);
            end
            %%Check Adducts
            %transform cell array to array
            rows=max(sum(cellfun(@numel,PossibleBaseMZindex)));
            idxMat=zeros(rows,length(mzVec));
            for k=1:length(mzVec)
                val=vertcat(PossibleBaseMZindex{:,k});
                idxMat(1:length(val),k)=val;
            end
            %find matching RT indices
            matchRT=cell(size(idxMat));
            for n=1:size(idxMat,2)
                for k=1:size(idxMat,1)
                    if idxMat(k,n)~=0
                        %check RT
                        [~,matchRT{k,n}]=ismembertol(RTs{n},RTs{idxMat(k,n)},1,'DataScale',1);
                    end
                end
            end
            %remove lists with only zeros
            check=cellfun(@(X) all(X(:)==0), matchRT);
            matchRT(check)={[]};
            idxMat(check)=0;
            %Preallocate storage variable of identified Isotopes variables
            AdductIndexCell = cellfun(@(X) false(size(X(:))), matchRT, 'UniformOutput', false);
            ColumnIndices = 1:1:size(matchRT,2);
            notemptyColumns = any(~cellfun(@isempty, matchRT),1);
            for n=ColumnIndices(notemptyColumns) %Isotope mz loop
                AdductEIC = EICMat(:,n);
                AdductRanges=ranges{1,n};
                RowIndices=1:1:size(matchRT,1);
                notEmptyRows = any(~cellfun(@isempty, matchRT(:,n)),2);
                for k=RowIndices(notEmptyRows) % possible Base mz loop
                    BaseEIC=EICMat(:,idxMat(k,n,1));    %base m/z EIC
                    BaseRanges=ranges{1,idxMat(k,n,1)}; %Base m/z peak ranges
                    CheckList=matchRT{k,n};             %List of Peaks that match Adduct rt
                    %preallocate variables
                    MatchIndices=1:1:size(CheckList,1);
                    notZeroMatch = CheckList ~=0;
                    for p=MatchIndices(notZeroMatch)   % loop over matching Peaks
                        RangeA=AdductRanges(p,:);
                        RangeB=BaseRanges(CheckList(p,1),:);
                        EICA=AdductEIC;
                        EICB=BaseEIC;
                        %cut vectors to size
                        minFullRange = min([RangeA;RangeB],[],"all");
                        maxFullRange = max([RangeA;RangeB],[],"all");
                        EICA = EICA(minFullRange:maxFullRange);
                        EICB = EICB(minFullRange:maxFullRange);
                        %check CosineSimilarity between adduct and MainPeak
                        AdductIndexCell{k,n}(p)=sum(EICB.*EICA)/(sqrt(sum(EICB.^2))*sqrt(sum(EICA.^2)))>=minCosSim ;
                    end
                end
            end
            %remove identified Isotopes
            [valuesFiltered,numRemoved]=obj.Removify(IntegrationResults,AdductIndexCell);
            obj.AdductFiltered = numRemoved;
            %remove Empty columns
            empt=cellfun(@isempty,valuesFiltered(2,:));
            valuesFiltered(:,empt)=[];
            obj.TempDataFileObj.ROIMat(:,empt) = [];
            obj.TempDataFileObj.ROImzVec(empt) = [];
        end

        function obj = filterIsotopes(obj)
            %% filter Isotopes removes isotope masses from centroided scans
            MASS_TOLERANCE = obj.withinFileMassTolerance;
            MASS_TOLERANCE_UNIT = obj.withinFileMassUnit;
            tempPeakData = obj.TempDataFileObj.ROICells;
            parfor iFile = 1:size(tempPeakData,1)
                tempPeakData{iFile,1} = scanIsotopeFilter(tempPeakData{iFile,1},MASS_TOLERANCE,MASS_TOLERANCE_UNIT);
            end
            obj.TempDataFileObj.ROICells = tempPeakData;
        end

        function obj = cutScansToSize(obj)
            %% removes scans outside specified time range and initializes TempDataFileObj
            START_TIME = obj.measurementStartTime;
            END_TIME = obj.measurementEndTime;
            tempPeakData = obj.RawDataFileObj.centroidedDataMS1;
            tempTimeData = obj.RawDataFileObj.timeDataMS1;

            parfor iFile = 1:size(tempPeakData,1)
                idx = tempTimeData{iFile,1} < START_TIME | tempTimeData{iFile,1} > END_TIME;
                tempPeakData{iFile,1}(idx)=[];
                tempTimeData{iFile,1}(idx)=[];
            end
            obj.TempDataFileObj.ROICells = tempPeakData;
            obj.TempDataFileObj.TimeCells = tempTimeData;
        end

        function obj = findRegionOfInterest(obj,modeFlag)
            %% AutoROI Performs fully automated ROI search and augmentation.
            switch modeFlag
                case "batch"
                    peakList = obj.TempDataFileObj.ROICells;
                    timeList = obj.TempDataFileObj.TimeCells;
                otherwise
                    peakList = obj.TempDataFileObj.ROICells(1,1);
                    timeList = obj.TempDataFileObj.TimeCells(1,1);
            end
            THRESHOLD = obj.roiThreshold;
            MIN_SIZE = obj.roiMinOccurrence;
            MASS_ERROR = obj.withinFileMassTolerance;
            MASS_ERROR_UNIT = obj.withinFileMassUnit;

            %preallocate cell arrays
            mzlist = cell(length(peakList),1);
            msRoiList = cell(length(peakList),1);

            %ROI search for every Sample
            parfor iFile = 1:length(peakList)
                currentPeak = peakList{iFile,1};
                currentTime = timeList{iFile,1};
                [mzlist{iFile,1},msRoiList{iFile,1},~] = ROIpeaks3(currentPeak,THRESHOLD,MASS_ERROR,MASS_ERROR_UNIT,MIN_SIZE,currentTime);
            end

            if numel(mzlist) > 1
                for iFile = 2:size(peakList,1)
                    [msRoiList{1,1},mzlist{1,1},timeList{1,1}] = MSroiaug3(msRoiList{1,1},msRoiList{iFile,1},mzlist{1,1},mzlist{iFile,1},MASS_ERROR,MASS_ERROR_UNIT,THRESHOLD,timeList{1,1},timeList{iFile,1});
                end
                msRoi_end = msRoiList{1,1};
                mzroi_end = mzlist{1,1};
                time_end = timeList{1,1};
            else %Skip Augmentation if only one Sample
                msRoi_end = msRoiList{1,1};
                mzroi_end = mzlist{1,1};
                time_end = timeList{1,1};
            end

            msRoi_end = msRoi_end-obj.roiThreshold; %subtract intensity threshold
            msRoi_end = max(msRoi_end,0); %set every negative intensity to 0

            if strcmp(modeFlag,"batch") %split and pad matrices
                outROI = mat2cell(msRoi_end,obj.nScan);
                outTime = mat2cell(time_end,obj.nScan);
                maxScan = max(obj.nScan);
                scanNumberArray = obj.nScan;
                parfor iFile = 1:size(outROI,1)
                    outROI{iFile,1} = padarray(outROI{iFile,1},maxScan-scanNumberArray(iFile,1),0,'post');
                    outTime{iFile,1} = padarray(outTime{iFile,1},maxScan-scanNumberArray(iFile,1),0,'post');
                end
            else
                outROI{1,1} = msRoi_end;
                outTime{1,1} = time_end;
            end
            obj.TempDataFileObj.ROICells = outROI;
            obj.TempDataFileObj.TimeCells = outTime;
            obj.TempDataFileObj.ROImzVec = mzroi_end;
        end

        function obj = alignMasses(obj,modeFlag)
            %% Applies mass alignment to centroided scans
            QUANTIL = obj.massAlignmentQuantil;
            ESTIMATION_METHOD = obj.massAlignmentEstimMethod;
            CORRECTION_METHOD = obj.massAlignmentCorrectionMethod;
            switch modeFlag
                case "preview"
                    peakCells = obj.TempDataFileObj.ROICells(1,1);
                otherwise
                    peakCells = obj.TempDataFileObj.ROICells;
            end
            parfor iFile = 1:size(peakCells,1)
                % perform Spectral Alignment
                [~, peakCells{iFile,1}] = mspalign(peakCells{iFile,1},'Quantile',QUANTIL,'EstimationMethod',ESTIMATION_METHOD,'CorrectionMethod',CORRECTION_METHOD,'ShowEstimation',false);
            end
            obj.TempDataFileObj.ROICells = peakCells;
        end

        function obj = correctBaseline(obj,modeFlag)
            %% applies baseline correction to ROIs
            switch modeFlag
                case "batch"
                    msRoi = obj.TempDataFileObj.ROICells;
                    time = obj.TempDataFileObj.TimeCells;
                otherwise
                    msRoi = obj.TempDataFileObj.ROICells(1,1);
                    time = obj.TempDataFileObj.TimeCells(1,1);
            end
            WINDOW_SIZE = obj.baselineWindowSize;
            STEP_SIZE = obj.baselineStepSize;
            REGRESSION = obj.baselineRegressionMethod;
            ESTIMATION = obj.baselineEstimationMethod;
            SMOOTHING = obj.baselineSmoothMethod;
            QUANTIL = obj.baselineQuantil;

            parfor iFile = 1:size(msRoi,1)
                oldSize = height(msRoi{iFile,1});
                %depad Array
                msRoiTemp = msRoi{iFile,1};
                [msRoiTemp,timeTemp] = depadArrays(msRoiTemp,time{iFile,1});
                msRoiTemp = msbackadj(timeTemp,msRoiTemp,'WindowSize',WINDOW_SIZE,'StepSize',STEP_SIZE,'RegressionMethod',REGRESSION,'EstimationMethod',ESTIMATION,'SmoothMethod',SMOOTHING,'QuantileValue',QUANTIL,'PreserveHeights',true);
                %remove negative, NaN and re-pad Array
                msRoiTemp = max(msRoiTemp,0);
                msRoiTemp(isnan(msRoiTemp)) = 0;
                [msRoi{iFile,1},time{iFile,1}] = repadArrays(msRoiTemp,timeTemp,oldSize);
                % set possible negative values to 0
                msRoi{iFile,1} = max(msRoi{iFile,1},0);
            end
            obj.TempDataFileObj.ROICells = msRoi;
            obj.TempDataFileObj.TimeCells = time;
        end

        function obj = smoothPeaks(obj,modeFlag)
            %% applies golay smoothing to ROIs
            switch modeFlag
                case "batch"
                    msRoi = obj.TempDataFileObj.ROICells;
                    time = obj.TempDataFileObj.TimeCells;
                otherwise
                    msRoi = obj.TempDataFileObj.ROICells(1,1);
                    time = obj.TempDataFileObj.TimeCells(1,1);
            end
            FRAME_SIZE = obj.smoothingFrameSize;
            DEGREE = obj.smoothingDegree;
            parfor iFile = 1:size(msRoi,1)
                %depad Array
                oldSize = height(msRoi{iFile,1});
                msRoiTemp = msRoi{iFile,1};
                [msRoiTemp,timeTemp] = depadArrays(msRoiTemp,time{iFile,1});
                msRoiTemp = mssgolay(timeTemp,msRoiTemp,'Span',FRAME_SIZE,'Degree',DEGREE);
                %remove negative, NaN and re-pad Array
                msRoiTemp = max(msRoiTemp,0);
                msRoiTemp(isnan(msRoiTemp))=0;
                [msRoi{iFile,1},time{iFile,1}] = repadArrays(msRoiTemp,timeTemp,oldSize);
                % set possible negative values to 0
                msRoi{iFile,1} = max(msRoi{iFile,1},0);
            end
            obj.TempDataFileObj.ROICells = msRoi;
            obj.TempDataFileObj.TimeCells = time;
        end

        function obj = alignPeaks(obj,modeFlag)
            %% applies peak alignment to ROIs
            switch modeFlag
                case "preview"
                    msRoi = obj.TempDataFileObj.ROICells(1,1);
                    time=obj.TempDataFileObj.TimeCells(1,1);
                case "batch"
                    msRoi = obj.TempDataFileObj.ROICells;
                    time=obj.TempDataFileObj.TimeCells;
            end
            WINDOW_SIZE_RATIO = obj.peakAlignmentWindowSizeRatio;
            ITERATION = obj.peakAlignmentIteration;
            GRID_STEPS = obj.peakAlignmentGridSteps;
            SEARCH_SPACE = obj.peakAlignmentSearchSpace;
            SHIFT_VALUES = [obj.peakAlignmentMaxShiftNegative*-1,obj.peakAlignmentMaxShiftPositive];
            PULSE_WIDTH = obj.peakAlignmentPulseWidth;
            maxScan = max(obj.nScan);

            % rearrange matrices
            [splitVar,~] = cellfun(@size,time);
            test= cellfun(@(x) sum(x~=0),time);
            [~,test] = max(test);
            timeVec = time{test};
            msRoi = vertcat(msRoi{:});
            time = vertcat(time{:});
            [~,id] = max(msRoi);
            referenceTimeArray = time(id);

            parfor iFile = 1:size(msRoi,2)
                currentRoi = reshape(msRoi(:,iFile),maxScan,[]);
                currentRoi = msalign(timeVec,currentRoi,referenceTimeArray(iFile),'MaxShift',SHIFT_VALUES,...
                    'WidthOfPulses',PULSE_WIDTH,'WindowSizeRatio',WINDOW_SIZE_RATIO,'Iterations',...
                    ITERATION,'GridSteps',GRID_STEPS,'SearchSpace',SEARCH_SPACE);
                currentRoi(isnan(currentRoi)) = 0; %remove possible NaN
                msTemp{1,iFile} = reshape(currentRoi,[],1);
            end
            msRoi = cell2mat(msTemp);
            obj.TempDataFileObj.ROICells = mat2cell(msRoi,splitVar);
            obj.TempDataFileObj.TimeCells = mat2cell(time,splitVar);
        end

        function obj = finalizeROI(obj)
            %% finalizes ROI output for integration
            % removes empty columns, applies padding and calculates scan
            % frequency
            time = obj.TempDataFileObj.TimeCells;
            msRoi = obj.TempDataFileObj.ROICells;
            paddedSize = zeros(size(msRoi));
            timeTemp = horzcat(time{:});
            timeTemp(any(timeTemp==0,2),:) = [];
            obj.scanFrequencySecond = mean(diff(timeTemp),'all');
            MAX_PEAK_WIDTH = round(obj.peakMaxWidth*1.5/obj.scanFrequencySecond);
            parfor iFile = 1:size(msRoi,1)
                msTemp = msRoi{iFile,1};
                msTemp = max(msTemp,0);
                msTemp(isnan(msTemp)) = 0;
                msTemp = padarray(msTemp,MAX_PEAK_WIDTH,0,'post');
                msRoi{iFile,1} = msTemp;
                time{iFile,1} = padarray(time{iFile,1},MAX_PEAK_WIDTH,0,'post');
                paddedSize(iFile) = length(time{iFile,1});
            end
            obj.nScanPadded = paddedSize;
            obj.nScan(length(paddedSize)+1:end) = [];
            msTemp = vertcat(msRoi{:});
            msTemp = max(msTemp,0);
            %remove empty columns
            id = all(msTemp == 0,1);
            msTemp(:,id) = [];
            if obj.useBlankSubtraction == true
                obj.TempDataFileObj.ROIMatBLK(:,id) = [];
            end
            obj.TempDataFileObj.ROImzVec(id) = [];
            obj.TempDataFileObj.ROIMat = sparse(msTemp);
            obj.TempDataFileObj.timeVec = round(vertcat(time{:}),1);
        end

        function IntegrationData = assignFileID(obj,IntegrationData)
            %% assignFileID
            % Assigns the original file identifier to each peak found by
            % integration
            MAX_FILE_SCAN = cumsum(obj.nScanPadded)';
            nFiles = length(obj.nScanPadded);

            parfor iFeature = 1:length(IntegrationData)
                fileID = zeros(size(IntegrationData(iFeature).peakLocation));
                peakLocation = IntegrationData(iFeature).peakLocation;
                for jFileID = 1:numel(fileID)
                    val = peakLocation(jFileID,1);
                    val = val < MAX_FILE_SCAN;
                    val = sum(val,2)-1;
                    fileID(jFileID,1) = abs(val-nFiles);
                end
                IntegrationData(iFeature).fileID = fileID;
            end
        end

        function OutputStruct = groupAndSampleScaling(obj,InputStruct)
            %% applies group and sample scaling to final output struct

            GROUP_SCALING_FACTORS = obj.groupScale;
            SAMPLE_SCALING_FACTORS = obj.sampleScale';

            OutputStruct = InputStruct;
            features = InputStruct.feature;

            if InputStruct.dataSize > 0
                parfor iFeature = 1:length(features)
                    %GroupScale
                    features(iFeature).peakHeights = features(iFeature).peakHeights/GROUP_SCALING_FACTORS;
                    features(iFeature).peakAreas = features(iFeature).peakAreas/GROUP_SCALING_FACTORS;
                    %SampleScale
                    features(iFeature).peakHeights = features(iFeature).peakHeights./SAMPLE_SCALING_FACTORS;
                    features(iFeature).peakAreas = features(iFeature).peakAreas./SAMPLE_SCALING_FACTORS;
                end
                OutputStruct.feature = features;
            end
        end

        function IntResults = filterPeaksFromIntegration(obj,IntResults,noiseArray)
            %% Filters identified peaks from CWT integration
            MinPWDataPoints = floor(obj.peakMinWidth/obj.scanFrequencySecond);
            MaxPWDataPoints = ceil(obj.peakMaxWidth/obj.scanFrequencySecond);
            maxSN = obj.minSignal2Noise;

            parfor iFeature = 1:length(IntResults)
                %check empty input
                if isempty(IntResults(iFeature).peakLocation)
                    continue
                end
                %% Peak filter
                %remove duplicate peaks
                uniquePeaks = unique([IntResults(iFeature).peakLocation,IntResults(iFeature).peakStartLocation,IntResults(iFeature).peakEndLocation,IntResults(iFeature).peakHeight,],'rows','stable');
                IntResults(iFeature).peakLocation = uniquePeaks(:,1);
                IntResults(iFeature).peakStartLocation = uniquePeaks(:,2);
                IntResults(iFeature).peakEndLocation = uniquePeaks(:,3);
                IntResults(iFeature).peakHeight = uniquePeaks(:,4);

                %preallocate index array
                idxToRemove = false(size(IntResults(iFeature).peakLocation));

                %remove peaks with wrong boundaries
                idToRemoveNew = IntResults(iFeature).peakStartLocation>=IntResults(iFeature).peakEndLocation;
                idxToRemove = idxToRemove | idToRemoveNew;

                %remove peaks with height = 0
                idToRemoveNew = IntResults(iFeature).peakHeight == 0;
                idxToRemove = idxToRemove | idToRemoveNew;

                %remove peaks with bad Peak asymmetry
                symmetry = (IntResults(iFeature).peakEndLocation - IntResults(iFeature).peakLocation)./(IntResults(iFeature).peakLocation - IntResults(iFeature).peakStartLocation);
                idToRemoveNew = symmetry<0.3 | symmetry>3;
                idxToRemove = idxToRemove | idToRemoveNew;

                %less than minimum peak width
                idToRemoveNew = IntResults(iFeature).peakEndLocation-IntResults(iFeature).peakStartLocation < MinPWDataPoints;
                IntResults(iFeature).minWidthFiltered=sum(idToRemoveNew);
                idxToRemove = idxToRemove | idToRemoveNew;

                %more than maximum peak width
                idToRemoveNew=IntResults(iFeature).peakEndLocation - IntResults(iFeature).peakStartLocation > MaxPWDataPoints;
                IntResults(iFeature).maxWidthFiltered=sum(idToRemoveNew);
                idxToRemove = idxToRemove | idToRemoveNew;

                %S/N peak rejection
                IntResults(iFeature).signal2Noise = IntResults(iFeature).peakHeight ./ noiseArray(iFeature);
                idToRemoveNew = IntResults(iFeature).signal2Noise < maxSN;
                IntResults(iFeature).signal2NoiseFiltered = sum(idToRemoveNew);
                idxToRemove = idxToRemove | idToRemoveNew;

                % remove identified peaks
                IntResults(iFeature).peakLocation(idxToRemove) = [];
                IntResults(iFeature).peakStartLocation(idxToRemove) = [];
                IntResults(iFeature).peakEndLocation(idxToRemove) = [];
                IntResults(iFeature).peakHeight(idxToRemove) = [];
                IntResults(iFeature).signal2Noise(idxToRemove) = [];
            end

            %calculate peak entropy and filter after first filter round
            % possible wrong peak bounderies causes errors
            IntResults = obj.calculatePeakEntropy(IntResults);

            %determine entropy bins
            if obj.useEntropyFilter == true
                allEntropy = vertcat(IntResults(:).entropy);
                [~,binedges] = histcounts(allEntropy,'BinMethod','auto');
                switch obj.entropyFilterStrength
                    case "lax"
                        medianEntropy = binedges(10);
                    case "medium"
                        medianEntropy = binedges(6);
                    case "strict"
                        medianEntropy = binedges(2);
                end

            else
                medianEntropy = 1;
            end

            parfor iFeature = 1:length(IntResults)
                %check empty input
                if isempty(IntResults(iFeature).peakLocation)
                    continue
                end
                %entropy peak rejection
                idToRemoveNew = IntResults(iFeature).entropy > medianEntropy;
                IntResults(iFeature).entropyFiltered = sum(idToRemoveNew);

                IntResults(iFeature).peakLocation(idToRemoveNew) = [];
                IntResults(iFeature).peakStartLocation(idToRemoveNew) = [];
                IntResults(iFeature).peakEndLocation(idToRemoveNew) = [];
                IntResults(iFeature).peakHeight(idToRemoveNew) = [];
                IntResults(iFeature).signal2Noise(idToRemoveNew) = [];
            end

            %remove features without peaks
            idToRemoveNew = false(length(IntResults),1);
            for iFeature = 1:length(IntResults)
                idToRemoveNew(iFeature) = isempty(IntResults(iFeature).peakLocation);
            end
            IntResults(idToRemoveNew) = [];

        end

        function OutputStruct = findOriginalMassScans(obj,InputStruct)
            %% gathers MS1 spectra for each feature and aligns them
            OutputStruct = InputStruct;
            allScans = obj.RawDataFileObj.profileDataMS1;
            % append all scans with spacers in between, to match processing
            % indices
            for iScan = 1:numel(obj.dataFile)
                temp = allScans{iScan,1};
                temp(obj.nScanPadded(iScan),1) = {[]};
                allScans{iScan,1} = temp;
            end
            allScans = vertcat(allScans{:});
            nFile = numel(obj.dataFile);

            parfor iScan = 1: length(InputStruct)
                spectra = cell(1,nFile);
                location = InputStruct(iScan).peakLocations;

                for jFile = 1:nFile
                    %check if borders contain NaN then skip iteration
                    if any(isnan(location(:,jFile)))
                        continue
                    end
                    %select spectra in peak range
                    scans = allScans(location(1,jFile));
                    %remove possible empty scans
                    scans(cellfun(@isempty, scans)) = [];
                    if ~isempty(scans)
                        %average scan
                        spectra{1,jFile} = alignSpectra(scans,"average","high","false");
                    else
                        spectra{1,jFile} = {[]};
                    end
                end
                OutputStruct(iScan).spectrumMS1 = spectra;
            end
        end

        function Output = occurrenceFilterFeatures(obj,Output)
            %% remove features with less peaks than required minimum from Output struct
            nFile = numel(obj.dataFile);
            minDataPoints = ceil(nFile*obj.minOccurence);
            FeatureStruct = Output.feature;
            numElements = zeros(length(FeatureStruct),1);

            for iFeature = 1:length(FeatureStruct)
                numElements(iFeature) = nnz(~isnan(FeatureStruct(iFeature).peakHeights));
            end
            idToRemove = numElements < minDataPoints;

            %sum number of removed peaks
            Output.occurenceFiltered = sum(numElements(idToRemove),"all");
            FeatureStruct(idToRemove) = [];

            Output.feature = FeatureStruct;
        end

        function Output = finalizeFeatureOutput(obj,Output)
            %% finalizes feature Output Struct
            % fills remaining fields: dataSize, featID
            % rounds retentionTime
            % not yet implemented: formula, adduct type, corrected mass

            Output.dataSize = length(Output.feature);
            %build average RetentionTime and featureID
            for iFeat = 1:Output.dataSize
                Output.feature(iFeat).retentionTime = mean(Output.feature(iFeat).retentionTimes,"all","omitmissing");
                Output.feature(iFeat).featID = Output.feature(iFeat).mass_measured + "Da@" + Output.feature(iFeat).retentionTime + "s_" + obj.groupName;
            end
            % separation Type specific tasks

            switch Output.separationType
                case "GC"
                    for iFeat = 1:Output.dataSize
                        % calculate formula
                        % GC -> database Search

                        % add Adduct type
                        Output.feature(iFeat).adductType = "M+";
                        % calculate corrected mass based on adduct Type or Database
                    end

                otherwise %ESI

                    % calculate formula
                    % ESI -> mass decomposition

                    % add Adduct type
                    % ESI check ms1 spectrum

                    % calculate corrected mass based on adduct Type or Database
            end

        end
        
        %% Data handling
        function obj = initializeStorageFile(obj)
            %check if file already exists
            if isfile(obj.RawDataFile)
                delete(obj.RawDataFile)
            end
            obj.RawDataFileObj = matfile(obj.RawDataFile,Writable=true);

            %predefine Variables in .mat file
            obj.RawDataFileObj.previewTICs = {[]};
            obj.RawDataFileObj.previewBPCs = {[]};
            obj.RawDataFileObj.previewTimes = {[]};

            obj.RawDataFileObj.profileDataMS1 = {[]};
            obj.RawDataFileObj.centroidedDataMS1 = {[]};
            obj.RawDataFileObj.timeDataMS1 = {[]};
            obj.RawDataFileObj.polarityMS1 = {[]};

            obj.RawDataFileObj.profileDataMS2 = {[]};
            obj.RawDataFileObj.centroidedDataMS2 = {[]};
            obj.RawDataFileObj.timeDataMS2 = {[]};
            obj.RawDataFileObj.polarityMS2 = {[]};
            obj.RawDataFileObj.precursorMass = {[]};
            obj.RawDataFileObj.molecularPrecursorMass = {[]};
            obj.RawDataFileObj.fragmentationEnergy = {[]};
            obj.RawDataFileObj.fragmentationType = {[]};
        end
    
        function obj = setOptimizationOptions(obj,optimizeMode,bayesOptions)

            switch optimizeMode

                case "MainOptions"
                    if ismember("intThresh",bayesOptions.Properties.VariableNames)
                        obj.roiThreshold = bayesOptions.intThresh;
                    end
                    if ismember("mzerror",bayesOptions.Properties.VariableNames)
                        obj.withinFileMassTolerance = bayesOptions.mzerror;
                    end
                    if ismember("minRoi",bayesOptions.Properties.VariableNames)
                        obj.roiMinOccurrence = bayesOptions.minRoi;
                    end
                    if ismember("minPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.minPeakWidth)
                        obj.peakMinWidth = bayesOptions.minPeakWidth;
                    end
                    if ismember("maxPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.maxPeakWidth)
                        obj.peakMaxWidth = bayesOptions.maxPeakWidth;
                    end
                    if ismember("minSN",bayesOptions.Properties.VariableNames)
                        obj.minSignal2Noise = bayesOptions.minSN;
                    end
                    if ismember("mzTol",bayesOptions.Properties.VariableNames)
                        obj.betweenFileMassTolerance = bayesOptions.mzTol;
                    end
                    if ismember("timeTol",bayesOptions.Properties.VariableNames)
                        obj.peakTimeTolerance = bayesOptions.timeTol;
                    end
                    if ismember("entropyFilter",bayesOptions.Properties.VariableNames)
                        obj.useEntropyFilter = bayesOptions.entropyFilter == "true";
                    end
                    if ismember("blankCorrection",bayesOptions.Properties.VariableNames)
                        obj.useBlankSubtraction = bayesOptions.blankCorrection == "true";
                    end
                    if ismember("contaminantFilter",bayesOptions.Properties.VariableNames)
                        obj.useContaminantFilter = bayesOptions.contaminantFilter == "true";
                    end
                    if ismember("isotopeFilter",bayesOptions.Properties.VariableNames)
                        obj.useIsotopeFilter = bayesOptions.isotopeFilter == "true";
                    end
                    if ismember("MSAlign",bayesOptions.Properties.VariableNames)
                        obj.useMassAlign = bayesOptions.MSAlign == "true";
                    end
                    if ismember("peakAlignment",bayesOptions.Properties.VariableNames)
                        obj.usePeakAlign = bayesOptions.peakAlignment == "true";
                    end
                    if ismember("baselineCorrection",bayesOptions.Properties.VariableNames)
                        obj.useBaselineCorrection = bayesOptions.baselineCorrection == "true";
                    end
                    if ismember("smoothing",bayesOptions.Properties.VariableNames)
                        obj.useSmoothing = bayesOptions.smoothing == "true";
                    end

                case "SubParameters"
                    if ismember("entropyStrength",bayesOptions.Properties.VariableNames)
                        obj.entropyFilterStrength = bayesOptions.entropyStrength;
                    end
                    %MSAlign parameters
                    if ismember("mzEstimMethod",bayesOptions.Properties.VariableNames)
                        obj.massAlignmentEstimMethod = bayesOptions.mzEstimMethod;
                    end
                    if ismember("mzCorrectionMethod",bayesOptions.Properties.VariableNames)
                        obj.massAlignmentCorrectionMethod = bayesOptions.mzCorrectionMethod;
                    end
                    if ismember("mzQuantil",bayesOptions.Properties.VariableNames)
                        obj.massAlignmentQuantil = bayesOptions.mzQuantil;
                    end
                    %PeakAlign parameters
                    if ismember("maxShiftneg",bayesOptions.Properties.VariableNames)
                        obj.peakAlignmentMaxShiftNegative = bayesOptions.maxShiftneg;
                    end
                    if ismember("maxShiftpos",bayesOptions.Properties.VariableNames)
                        obj.peakAlignmentMaxShiftPositive = bayesOptions.maxShiftpos;
                    end
                    if ismember("pulseWidth",bayesOptions.Properties.VariableNames)
                        obj.peakAlignmentPulseWidth = bayesOptions.pulseWidth;
                    end
                    if ismember("iterations",bayesOptions.Properties.VariableNames)
                        obj.peakAlignmentIteration = bayesOptions.iterations;
                    end
                    if ismember("searchSpace",bayesOptions.Properties.VariableNames)
                        obj.peakAlignmentSearchSpace = bayesOptions.searchSpace;
                    end
                    if ismember("gridSteps",bayesOptions.Properties.VariableNames)
                        obj.peakAlignmentGridSteps = bayesOptions.gridSteps;
                    end
                    %baseline parameters
                    %always set smoothing to none
                    obj.baselineSmoothMethod = "none";
                    if ismember("windowSize",bayesOptions.Properties.VariableNames)
                        obj.baselineWindowSize = bayesOptions.windowSize;
                    end
                    if ismember("stepSize",bayesOptions.Properties.VariableNames)
                        obj.baselineStepSize = bayesOptions.stepSize;
                    end
                    if ismember("regressionMethod",bayesOptions.Properties.VariableNames)
                        obj.baselineRegressionMethod = bayesOptions.regressionMethod;
                    end
                    if ismember("estimationMethod",bayesOptions.Properties.VariableNames)
                        obj.baselineEstimationMethod = bayesOptions.estimationMethod;
                    end
                    if ismember("quantile",bayesOptions.Properties.VariableNames)
                        obj.baselineQuantil = bayesOptions.quantile;
                    end
                    %Smoothing parameters
                    if ismember("frameSize",bayesOptions.Properties.VariableNames)
                        obj.smoothingFrameSize = bayesOptions.frameSize;
                    end
                    if ismember("polyDegree",bayesOptions.Properties.VariableNames)
                        obj.smoothingDegree = bayesOptions.polyDegree;
                    end


                case {"Full","Custom"}
                    if ismember("intThresh",bayesOptions.Properties.VariableNames)
                        obj.roiThreshold = bayesOptions.intThresh;
                    end
                    if ismember("mzerror",bayesOptions.Properties.VariableNames)
                        obj.withinFileMassTolerance = bayesOptions.mzerror;
                    end
                    if ismember("minRoi",bayesOptions.Properties.VariableNames)
                        obj.roiMinOccurrence = bayesOptions.minRoi;
                    end
                    if ismember("minPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.minPeakWidth)
                        obj.peakMinWidth = bayesOptions.minPeakWidth;
                    end
                    if ismember("maxPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.maxPeakWidth)
                        obj.peakMaxWidth = bayesOptions.maxPeakWidth;
                    end
                    if ismember("minSN",bayesOptions.Properties.VariableNames)
                        obj.minSignal2Noise = bayesOptions.minSN;
                    end
                    if ismember("mzTol",bayesOptions.Properties.VariableNames)
                        obj.betweenFileMassTolerance = bayesOptions.mzTol;
                    end
                    if ismember("timeTol",bayesOptions.Properties.VariableNames)
                        obj.peakTimeTolerance = bayesOptions.timeTol;
                    end
                    if ismember("entropyFilter",bayesOptions.Properties.VariableNames)
                        obj.useEntropyFilter = bayesOptions.entropyFilter == "true";
                        if bayesOptions.entropyFilter == "true"
                            if ismember("entropyStrength",bayesOptions.Properties.VariableNames)
                                obj.entropyFilterStrength = bayesOptions.entropyStrength;
                            end
                        end
                    end
                    if ismember("blankCorrection",bayesOptions.Properties.VariableNames)
                        obj.useBlankSubtraction = bayesOptions.blankCorrection == "true";
                    end
                    if ismember("contaminantFilter",bayesOptions.Properties.VariableNames)
                        obj.useContaminantFilter = bayesOptions.contaminantFilter == "true";
                    end
                    if ismember("isotopeFilter",bayesOptions.Properties.VariableNames)
                        obj.useIsotopeFilter = bayesOptions.isotopeFilter == "true";
                    end
                    if ismember("MSAlign",bayesOptions.Properties.VariableNames)
                        obj.useMassAlign = bayesOptions.MSAlign == "true";
                    end
                    if ismember("peakAlignment",bayesOptions.Properties.VariableNames)
                        obj.usePeakAlign = bayesOptions.peakAlignment == "true";
                    end
                    if ismember("baselineCorrection",bayesOptions.Properties.VariableNames)
                        obj.useBaselineCorrection = bayesOptions.baselineCorrection == "true";
                    end
                    if ismember("smoothing",bayesOptions.Properties.VariableNames)
                        obj.useSmoothing = bayesOptions.smoothing == "true";
                    end
                    %MS Alignment
                    if ismember("MSAlign",bayesOptions.Properties.VariableNames)
                        obj.useMassAlign = bayesOptions.MSAlign == "true";
                        if bayesOptions.MSAlign == "true"
                            if ismember("mzEstimMethod",bayesOptions.Properties.VariableNames)
                                obj.massAlignmentEstimMethod = bayesOptions.mzEstimMethod;
                            end
                            if ismember("mzCorrectionMethod",bayesOptions.Properties.VariableNames)
                                obj.massAlignmentCorrectionMethod = bayesOptions.mzCorrectionMethod;
                            end
                            if ismember("mzQuantil",bayesOptions.Properties.VariableNames)
                                obj.massAlignmentQuantil = bayesOptions.mzQuantil;
                            end
                        end
                    end
                    if ismember("peakAlignment",bayesOptions.Properties.VariableNames)
                        obj.usePeakAlign = bayesOptions.peakAlignment == "true";
                        if bayesOptions.peakAlignment == "true"
                            %Peak Alignment
                            if ismember("maxShiftneg",bayesOptions.Properties.VariableNames)
                                obj.peakAlignmentMaxShiftNegative = bayesOptions.maxShiftneg;
                            end
                            if ismember("maxShiftpos",bayesOptions.Properties.VariableNames)
                                obj.peakAlignmentMaxShiftPositive = bayesOptions.maxShiftpos;
                            end
                            if ismember("pulseWidth",bayesOptions.Properties.VariableNames)
                                obj.peakAlignmentPulseWidth = bayesOptions.pulseWidth;
                            end
                            if ismember("iterations",bayesOptions.Properties.VariableNames)
                                obj.peakAlignmentIteration = bayesOptions.iterations;
                            end
                            if ismember("searchSpace",bayesOptions.Properties.VariableNames)
                                obj.peakAlignmentSearchSpace = bayesOptions.searchSpace;
                            end
                            if ismember("gridSteps",bayesOptions.Properties.VariableNames)
                                obj.peakAlignmentGridSteps = bayesOptions.gridSteps;
                            end
                        end
                    end
                    %Baseline Correction
                    if ismember("baselineCorrection",bayesOptions.Properties.VariableNames)
                        obj.useBaselineCorrection = bayesOptions.baselineCorrection == "true";
                        if bayesOptions.baselineCorrection == "true"
                            %always set smoothing to none
                            obj.baselineSmoothMethod = "none";
                            if ismember("windowSize",bayesOptions.Properties.VariableNames)
                                obj.baselineWindowSize = bayesOptions.windowSize;
                            end
                            if ismember("stepSize",bayesOptions.Properties.VariableNames)
                                obj.baselineStepSize = bayesOptions.stepSize;
                            end
                            if ismember("regressionMethod",bayesOptions.Properties.VariableNames)
                                obj.baselineRegressionMethod = bayesOptions.regressionMethod;
                            end
                            if ismember("estimationMethod",bayesOptions.Properties.VariableNames)
                                obj.baselineEstimationMethod = bayesOptions.estimationMethod;
                            end
                            if ismember("smoothingMethod",bayesOptions.Properties.VariableNames)
                                obj.baselineSmoothMethod = bayesOptions.smoothingMethod;
                            end
                            if ismember("quantile",bayesOptions.Properties.VariableNames)
                                obj.baselineQuantil = bayesOptions.quantile;
                            end
                        end
                    end
                    % Smoothing
                    if ismember("smoothing",bayesOptions.Properties.VariableNames)
                        obj.useSmoothing = bayesOptions.smoothing == "true";
                        if bayesOptions.smoothing == "true"
                            if ismember("frameSize",bayesOptions.Properties.VariableNames)
                                obj.smoothingFrameSize = bayesOptions.frameSize;
                            end
                            if ismember("polyDegree",bayesOptions.Properties.VariableNames)
                                obj.smoothingDegree = bayesOptions.polyDegree;
                            end
                        end
                    end
            end
        end
    
    end

    methods (Static)

        function IntegrationStruct = calculatePeakEntropy(IntegrationStruct)
            % Calculates Peak entropy for all peaks
            for iFeature = 1:length(IntegrationStruct)
                %check for no peaks, then skip iteration
                if ~isempty(IntegrationStruct(iFeature).peakStartLocation)
                    derivativ = diff(IntegrationStruct(iFeature).XIC(:,2));
                    probability = zeros(size(IntegrationStruct(iFeature).peakLocation));
                    parfor iPeak = 1:numel(probability)
                        %extract peak range
                        Peak = derivativ(IntegrationStruct(iFeature).peakStartLocation(iPeak,:):IntegrationStruct(iFeature).peakEndLocation(iPeak,:));
                        maxidx = IntegrationStruct(iFeature).peakLocation(iPeak)-IntegrationStruct(iFeature).peakStartLocation(iPeak);
                        % check normal or variant point , variant point = 1
                        premax = Peak(1:maxidx-1) < 0;
                        postmax = Peak(maxidx+1:end) > 0;
                        VarPoints = [premax; false; postmax];
                        %calculate probability of variant point
                        probability(iPeak,1) = sum(VarPoints)/numel(VarPoints);
                    end
                    %calculate entropy
                    peakEntropy = -probability.*log2(probability)-(1-probability).*log2(1-probability);
                    peakEntropy(isnan(peakEntropy)) = 0;
                    %store values
                    IntegrationStruct(iFeature).entropy = peakEntropy;
                else
                    continue
                end
            end
        end

        function OutArray = fileSortPeaks(InArray)
            %% sort struct contents to original file

            %check if GC or LC/CE
            isGC = isscalar(InArray);

            %preallocate output Feature struct
            OutArray = InArray;

            nFiles = max(vertcat(InArray(:).fileID));
            %sort
            parfor iFeature = 1:length(InArray)
                %preallocate temporary cells
                massCell = cell(1,nFiles);
                peakLcell = cell(1,nFiles);
                peakRTcell = cell(1,nFiles);
                peakSLcell = cell(1,nFiles);
                peakELcell = cell(1,nFiles);
                peakHcell = cell(1,nFiles);
                peakAcell = cell(1,nFiles);
                entroCell = cell(1,nFiles);
                s2ncell = cell(1,nFiles);
                sMS2cell = cell(1,nFiles);
                fileidcell = cell(1,nFiles);

                id = InArray(iFeature).fileID;
                for jFile = 1:nFiles
                    % sort peaks into respective file cells
                    idx = id == jFile;
                    peakLcell{jFile} = InArray(iFeature).peakLocation(idx);
                    peakRTcell{jFile} = InArray(iFeature).peakRetentionTime(idx);
                    peakSLcell{jFile} = InArray(iFeature).peakStartLocation(idx);
                    peakELcell{jFile} = InArray(iFeature).peakEndLocation(idx);
                    peakHcell{jFile} = InArray(iFeature).peakHeight(idx);
                    peakAcell{jFile} = InArray(iFeature).peakArea(idx);
                    entroCell{jFile} = InArray(iFeature).entropy(idx);
                    s2ncell{jFile} = InArray(iFeature).signal2Noise(idx);
                    fileidcell{jFile} = InArray(iFeature).fileID(idx);
                    if isGC
                        massCell{jFile} = InArray(iFeature).mass(idx);
                        sMS2cell{jFile} = InArray(iFeature).spectrumMS2(idx);
                    end
                end
                %store temp cells into output
                OutArray(iFeature).peakLocation = peakLcell;
                OutArray(iFeature).peakRetentionTime = peakRTcell;
                OutArray(iFeature).peakStartLocation = peakSLcell;
                OutArray(iFeature).peakEndLocation = peakELcell;
                OutArray(iFeature).peakHeight = peakHcell;
                OutArray(iFeature).peakArea = peakAcell;
                OutArray(iFeature).entropy = entroCell;
                OutArray(iFeature).signal2Noise = s2ncell;
                OutArray(iFeature).fileID = fileidcell;

                if isGC
                    OutArray(iFeature).mass = massCell;
                    OutArray(iFeature).spectrumMS2 = sMS2cell;
                end
            end
        end

        function integrationStruct = finalizeIntegrationOutput(integrationStruct,times)
            %% Performs Integration of found Peaks and gathers retention times
            % get peak area and final retention time

            %check for empty struct
            if isempty(integrationStruct)
                return
            end

            parfor iFeature = 1:length(integrationStruct)
                eic = full(integrationStruct(iFeature).XIC);
                areas = zeros(size(integrationStruct(iFeature).peakLocation));
                for jPeak = 1:numel(areas)
                    areas(jPeak,1) = trapz(eic(integrationStruct(iFeature).peakStartLocation(jPeak,1):integrationStruct(iFeature).peakEndLocation(jPeak,1)));
                end
                integrationStruct(iFeature).peakArea = areas;
                integrationStruct(iFeature).peakRetentionTime = times(integrationStruct(iFeature).peakLocation);
            end
        end

        function Output = confirmSameFeatureByMS2(Output)
            %% calculates the composit score within each feature
            % features with scores < 650 are split into a new feature

            FeatureStruct = Output.feature;
            nFeat = length(FeatureStruct);
            
            MINIMUM_SCORE = 650;

            for iFeat = 1:nFeat
                spectra = FeatureStruct(iFeat).spectrumMS2;
                if sum(~cellfun("isempty",spectra)) <= 1 %only one file with peak or no spectra
                    continue
                end

                %build index to original file
                originalFileID = 1:numel(spectra);
                originalFileID(cellfun(@isempty,spectra)) = [];
                spectra(cellfun(@isempty,spectra)) = [];

                % clean scans
                for jSpectrum = 1:width(spectra)
                    data = spectra{1,jSpectrum}{:};
                    if isempty(data)
                        continue
                    end
                    idx = data(:,2) < 0.05;
                    data(idx,:) = [];
                    spectra{1,jSpectrum} = data;
                end

                alingedSpectra = alignSpectra(spectra,"normal","low","true");
                compoundScores = scoresWithinSet(alingedSpectra);
                % rebuild original file
                for jFile = 1:numel(originalFileID)
                    compoundScores(compoundScores==jFile) = originalFileID(jFile);
                end

                if all(compoundScores(:,1) >= 650)
                    continue
                else %split feature
                    idtoKeep = compoundScores(compoundScores(:,1) >= MINIMUM_SCORE,2:3);
                    idtoKeep = unique(idtoKeep);
                    idtoSplit = compoundScores(compoundScores(:,1) < MINIMUM_SCORE,2:3);
                    idtoSplit = unique(idtoSplit);
                    id = any(idtoSplit == idtoKeep,1);
                    idtoSplit(id) = [];

                    %check number of Peaks to remove
                    if isempty(idtoSplit)       %none, because of overlap
                        continue
                    elseif isempty(idtoKeep)    %all, keep first entry remove the rest
                        idtoSplit(1) = [];
                    end

                    splitFeatures = SplitFeature(FeatureStruct(iFeat),idtoSplit);
                    %append
                    FeatureStruct(iFeat) = splitFeatures(1);
                    FeatureStruct = [FeatureStruct;splitFeatures(2:end)];
                end
            end
            Output.feature = FeatureStruct;
        end

        function Output = confirmSameFeatureByIsotopeDistribution(Output)
            %% WIP
            featureStruct = Output.feature;
            numFeats = length(featureStruct);

            for n = 1:numFeats
                spectra = featureStruct(n).spectrumMS1;
                if sum(~cellfun("isempty",spectra)) <= 1 %only one file with peak or no spectra
                    continue
                end

                %build index to original file
                originalFileID = 1:numel(spectra);
                originalFileID(cellfun(@isempty,spectra)) = [];
                spectra(cellfun(@isempty,spectra)) = [];

                % clean scans
                for j = 1:width(spectra)
                    data = spectra{1,j};
                    if isempty(data)
                        continue
                    end
                    idx = data(:,2) < 0.05;
                    data(idx,:) = [];
                    spectra{1,j} = data;
                end

                alingedSpectra = AlignSpectra(spectra,"average","high","true");
                CompoundScores = ScoresWithinSet(alingedSpectra);
                % rebuild original file
                for file = 1:numel(originalFileID)
                    CompoundScores(CompoundScores==file) = originalFileID(file);
                end

                if all(CompoundScores(:,1) >= 650)
                    continue
                else %split feature

                    idtoKeep = CompoundScores(CompoundScores(:,1) >= 650,2:3);
                    idtoKeep = unique(idtoKeep);
                    idtoSplit = CompoundScores(CompoundScores(:,1) < 650,2:3);
                    idtoSplit = unique(idtoSplit);
                    id = any(idtoSplit == idtoKeep,1);
                    idtoSplit(id) = [];

                    %check number of Peaks to remove
                    if isempty(idtoSplit)       %none, because of overlap
                        continue
                    elseif isempty(idtoKeep)    %all, keep first entry remove the rest
                        idtoSplit(1) = [];
                    end

                    splitFeatures = SplitFeature(featureStruct(n),idtoSplit);
                    %append
                    featureStruct(n) = splitFeatures(1);
                    featureStruct = [featureStruct;splitFeatures(2:end)];
                end
            end
            Output.feature = featureStruct;
        end

    end
end