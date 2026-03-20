classdef RawData
    %% RawData superclass object for AriumMS file processing
    % Processing includes:
    % - Loading of MS data files
    % - Pre processing
    % - ROI search
    % - Post processing
    % - Feature extraction
    properties
        %% Processing Parameters
        fileName    (:,1) string
        dataFile    (:,1) string
        fileType    (:,1) string
        groupName   (1,1) string
        ionisationType (1,1) string {mustBeMember(ionisationType,["Soft","Hard"])} = "Soft"
        % Main Processing Options
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
        baselineQuantil             (1,1) double {mustBeInRange(baselineQuantil,0,1)} = 0.1
        % Golay Parameters
        smoothingFrameSize  (1,1) double {mustBeInteger,mustBePositive} = 20
        smoothingDegree     (1,1) double {mustBeInteger,mustBePositive} = 2
        % Internal Standard Data
        nInternalStandard     (1,1) double {mustBeInteger,mustBePositive} = 1
        internalStandardData  (:,3) double
        useISMassCorrection   (1,1) logical = false
        % Adduct Parameters
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
        RawDataFileObj      (1,:)
        TempDataFile        string
        TempDataFileObj     (1,1)
        ROIDataFile         string
        ROIDataFileObj      (1,1)
        % processing variables
        nScan              (:,1) double {mustBeInteger,mustBePositive}
        nScanPadded        (:,1) double {mustBeInteger,mustBePositive}
        occurenceFiltered  (1,1) double {mustBeInteger} = 0
        isSampleFile       (:,1) logical = true
        %FileInfos
        DataInfo             (1,:) struct
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
            else
                obj.mainWindow = uifigure;
            end
            obj.groupName = groupName;
            obj.ROIDataFile = tempname +".mat";

        end
        %% Handling MS files

        function obj = readData(obj)
            files = obj.dataFile;
            nFile = size(files,1);
            obj.isSampleFile = true(nFile,1);
            %preallocation
            MSData = struct('file',[],...
                'startTimeStamp',[],...
                'nSpectra',[],...
                'spectraMS1',[],...
                'spectraMS2',[]);
            
             FileInfo =  struct('numberOfScansMS1',[],...
                'numberOfScansMSn',[],...
                'startTime',[],...
                'endTime',[],...
                'scanFrequencyMS1',[],...
                'scanFrequencyMSn',[]);

            parfor iFile = 1:nFile
                % read file into struct
                currentFile = importMZML(files{iFile});

                % gather scan infos
                FileInfo(iFile).numberOfScansMS1 = numel([currentFile.spectraMS1]);
                FileInfo(iFile).numberOfScansMSn = numel([currentFile.spectraMS2]);
                FileInfo(iFile).startTime = currentFile.spectraMS1(1).rt;
                FileInfo(iFile).endTime = currentFile.spectraMS1(end).rt;
                
                MSData(iFile) = currentFile;
            end

            % calculate Scan Frequency [Hz]
            scanFrequencyHertz = [FileInfo.numberOfScansMS1]./([FileInfo.endTime]-[FileInfo.startTime]);
            scanFrequencyHertz = num2cell(scanFrequencyHertz);
            [FileInfo.scanFrequencyMS1] = scanFrequencyHertz{:};
            scanFrequencyHertz = [FileInfo.numberOfScansMSn]./([FileInfo.endTime]-[FileInfo.startTime]);
            scanFrequencyHertz = num2cell(scanFrequencyHertz);
            [FileInfo.scanFrequencyMSn] = scanFrequencyHertz{:};
            obj.DataInfo = FileInfo;
            obj.measurementStartTime = round(min([FileInfo.startTime]),1);
            obj.measurementEndTime = round(max([FileInfo.endTime]),1);

            % sort MS data based on start time stamp
            [~, idx] = sort([MSData.startTimeStamp]);    % ascending
            % apply sort to struct
            MSData = MSData(idx);

            % %% Change timestamps for test sequence data
            % MSData(4).startTimeStamp(1) = '30-Aug-2023 20:00:00';
            % MSData(5).startTimeStamp(1) = '30-Aug-2023 23:50:01';
            % MSData(6).startTimeStamp(1) = '31-Aug-2023 03:00:52';
            % [~, idx] = sort([MSData.startTimeStamp]);    % ascending
            % % apply sort to struct
            % MSData = MSData(idx);
            % obj.fileName = {MSData.file}';

            % determine file type (blank,QC,sample)
            obj = obj.determineFileType;
            
            obj.RawDataFileObj = MSData;
        end

        %% Data Processing
        function [Output,obj] = extractFeaturesFromMassData(obj,varargin)
            %check if old results exist and delete them
            if isfile(obj.ROIDataFile)
                delete(obj.ROIDataFile)
            end

            if numel(varargin) == 2
                mode = varargin{1};
                bayesOptions = varargin{2};
                obj = obj.setOptimizationOptions(mode,bayesOptions);
            elseif isscalar(varargin)
                error("Wrong number of inputs")
            else
                title = "Processing " + obj.groupName;
                progressBar = uiprogressdlg(obj.mainWindow,"Title",title,"Message","Preparation",Value=0);
            end
            
            obj = obj.initializeTemporaryFile;
            
            %remove scans outside RT range
            progressBar.Value = 0.1;
            progressBar.Message = "Loading data";
            obj = obj.cutScansToSize;
            obj.nScan = cellfun(@numel,obj.TempDataFileObj.TimeCells);
            progressBar.Value = 0.33;
            progressBar.Message = "Processing";
            switch obj.ionisationType
                case "Hard"
                    [Output,obj] = obj.extractFeatures_HardIonisation(progressBar);
                case "Soft"
                    [Output,obj] = obj.extractFeatures_SoftIonisation(progressBar);
            end

            %processing cleanup
            if ~exist("mode","var") %save ROI results if batch mode
                %build ResultDataFile
                obj.ROIDataFileObj = matfile(obj.ROIDataFile,Writable=true);
                %store Results
                obj.ROIDataFileObj.ROIMat = obj.TempDataFileObj.ROIMat;
                obj.ROIDataFileObj.ROIMatBLK  = obj.TempDataFileObj.ROIMatBLK;
                obj.ROIDataFileObj.ROImzVec = obj.TempDataFileObj.ROImzVec;
                obj.ROIDataFileObj.timeVec  = obj.TempDataFileObj.timeVec;
            end
            %delete Temporary file
            delete(obj.TempDataFile)
            obj.TempDataFile = "";
            close(progressBar)
        end
        %% Soft Ionisation Functions
        function [Output,obj] = extractFeatures_SoftIonisation(obj,progressBar)
            % remove isotopes
            if obj.useIsotopeFilter == true
                progressBar.Message = "Removing Isotopes";
                obj = obj.filterIsotopes;
                progressBar.Value = progressBar.Value + 0.05;
            end

            if obj.useMassAlign == true
                progressBar.Message = "Aligning MS Scans";
                obj = obj.alignMasses("batch");
                progressBar.Value = 0.4;
            end

            % ROI Search
            progressBar.Message = "Searching for ROIs";
            obj = obj.findRegionOfInterest("batch");
            progressBar.Value = 0.5;

            %Common Contaminant filter
            if obj.useContaminantFilter == true
                progressBar.Message = "Removing Contaminants";
                obj = obj.removeContaminants;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Baseline Correction
            if obj.useBaselineCorrection == true
                progressBar.Message = "Correcting Baseline";
                obj = obj.correctBaseline("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Smoothing
            if obj.useSmoothing == true
                progressBar.Message = "Smoothing Peaks";
                obj = obj.smoothPeaks("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Peak Align
            if obj.usePeakAlign == true && nData > 1
                progressBar.Message = "Aligning Peaks";
                obj = obj.alignPeaks("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end
            
            % Blank Correction
            if any(contains(obj.fileType,"blank")) % Separate Blank data from Sample data and subtract
                progressBar.Message = "Subtracting Blank";
                obj = obj.blankCorrection;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % pad arrays with Maximum peak width*1.5 Scans to eliminate
            % integration interference between matrices
            obj = obj.finalizeROI;

            clearvars -except obj progressBar
            %% Integration Stage
            % Find and Integrate IS separate
            if obj.useInternalStandard == true
                progressBar.Message = "Searching for Internal Standard";
                obj = obj.identifyInternalStandard;
                if ~isempty(obj.internalStandardData)
                    obj = obj.internalStandardNormalization;
                end
                progressBar.Value = progressBar.Value + 0.05;
            end

            % mass correction
            if obj.useISMassCorrection == true && ~isempty(obj.internalStandardMassDelta)
                progressBar.Message = "Performing IS mass correction";
                obj = obj.massCorrectionByInternalStandard;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Integrate all Peaks
            IDX = true(1,size(obj.TempDataFileObj.ROIMat,2));
            progressBar.Message = "Integrating Peaks";
            IntegrationData = obj.findPeaks_Soft(IDX);
            progressBar.Value = 0.9;

            progressBar.Message = "Processing found Features";
            IntegrationData = obj.assignFileID(IntegrationData);
            IntegrationData = obj.fileSortPeaks(IntegrationData);


            %%%%%%%
            % % remove adducts
            % if obj.useAdductFilter == true
            %     [IntegrationData,obj] = obj.filterAdducts(IntegrationData);
            % end
            %%%%%%%

            % build output Struct and store number of removed peaks per
            % filter
            Output = obj.initializeOutputStruct(IntegrationData);

            % Build Storage Arrays and filter by number of occurrences
            switch obj.ionisationType
                case "Hard"
                    [FeatureData,obj] = obj.buildFeatureArray_Hard(IntegrationData);
                case "Soft"
                    [FeatureData,obj] = obj.buildFeatureArray_Soft(IntegrationData);
            end
           
            clearvars IntegrationData

            % apply scaling
            if obj.useScaling == true
                progressBar.Message = "Apply scaling";
                FeatureData = obj.groupAndSampleScaling(FeatureData);
            end

            %fill remaining fields
            Output = obj.finalizeOutputStruct(Output,FeatureData);

            progressBar.Message = "Group processing successful";
            progressBar.Value = 1;

            obj.Output = Output;
            delete(progressBar)
        end

        function obj = determineFileType(obj)
            % determine file type (blank, QC sample) from file name
            dataType = strings(size(obj.fileName));
            dataType(1:end) = "sample";

            name = obj.fileName;

            % check filename for 'Blank' or 'QC' (case-insensitive)
            idxBlankProcess = contains(name,["ProcessBlank","Blank"],"IgnoreCase",true);
            idxBlankSystem = contains(name,["BW","SystemBlank"],"IgnoreCase",true);
            idxQC = contains(name,["50uM_MS","QC","Control","PooledQC"],"IgnoreCase",true);
            idxSystemSuitability = contains(name,"Suitability","IgnoreCase",true);

            dataType(idxQC) = "QC";
            dataType(idxSystemSuitability) = "suitability QC";
            dataType(idxBlankProcess) = "process blank";
            dataType(idxBlankSystem) = "suitability blank";
            obj.fileType = dataType;

        end

        function IntegrationResults = findPeaks_Soft(obj,Index)
            %gather data
            roiMat = obj.TempDataFileObj.ROIMat;
            roiMat = roiMat(:,Index);
            mzValueArray = obj.TempDataFileObj.ROImzVec;
            timeArray = obj.TempDataFileObj.timeVec;
            minPeakWidthScans = floor(obj.peakMinWidth/obj.scanFrequencySecond);

            % calculate EIC derivatives and store as sparse
            smoothed = smoothdata(roiMat,"gaussian","omitnan","SmoothingFactor",0.1);
            Diff2 = zeros(length(timeArray),size(roiMat,2));
            Diff2(1:end-2,:) = diff(smoothed,2);
            nMass = size(roiMat,2);

            % prepare wavelet filter-bank
            FilterBank = cwtfilterbank("SignalLength",size(Diff2,1), ...
                "WaveletParameters",[3 4], ...
                "VoicesPerOctave",8, ...
                "SamplingPeriod",seconds(obj.scanFrequencySecond), ...
                "PeriodLimits",[seconds(obj.peakMinWidth) seconds(obj.peakMaxWidth)]);

            %preallocate storage struct
            IntegrationResults = struct( ...
                "mass",[], ...
                "peakLocation",[], ...
                "peakRetentionTime",[], ...
                "peakStartLocation",[], ...
                "peakEndLocation",[], ...
                "peakHeight",[], ...
                "peakArea",[], ...
                "entropy",[], ...
                "signal2Noise",[], ...
                "minWidthFiltered",[], ...
                "maxWidthFiltered",[], ...
                "entropyFiltered",[], ...
                "signal2NoiseFiltered",[], ...
                "XIC",[], ...
                "fileID",[]);

            IntegrationResults = repmat(IntegrationResults,nMass,1);

            parfor iMass = 1:nMass
                peakData = continuosWaveletPeakPicking(Diff2(:,iMass),smoothed(:,iMass),FilterBank);
                currentEIC = roiMat(:,iMass);
                % Correct Peak Borders
                peakData = correctPeakData(peakData,currentEIC,smoothed(:,iMass),minPeakWidthScans);
                IntegrationResults(iMass).mass = mzValueArray(iMass);
                IntegrationResults(iMass).peakLocation = peakData(:,1);
                IntegrationResults(iMass).peakStartLocation = peakData(:,2);
                IntegrationResults(iMass).peakEndLocation = peakData(:,3);
                IntegrationResults(iMass).peakHeight = peakData(:,4);
                %store EIC
                IntegrationResults(iMass).XIC = [timeArray,currentEIC];
            end
            %filter found peaks
            noise = std(roiMat-smoothed);
            IntegrationResults = obj.filterPeaksFromIntegration(IntegrationResults,noise);
            IntegrationResults = obj.integratePeaks(IntegrationResults,timeArray);
        end

        function [FeatureResults,obj] = buildFeatureArray_Soft(obj,IntegrationResults,varargin)

            FeatureResults = rmfield(IntegrationResults,["minWidthFiltered","maxWidthFiltered","entropyFiltered","signal2NoiseFiltered"]);
            
            if isscalar(varargin)
                isISIntegration = true;
            else
                isISIntegration = false;
            end
            %check for empty IntegrationResults
            if ~isempty(IntegrationResults)

                nFiles = numel(obj.RawDataFileObj);

                EmptyStruct = struct(...
                    "ID",strings,...
                    "mz_measured",[],...
                    "rt",NaN,...
                    "Int",NaN(1,nFiles),...
                    "Area",NaN(1,nFiles),...
                    "mz_s",NaN(1,nFiles),...
                    "mz_e",NaN(1,nFiles),...
                    "rt_s",NaN(1,nFiles),...
                    "rt_e",NaN(1,nFiles),...
                    "PWMD",NaN(1,nFiles),...
                    "PWTD",NaN(1,nFiles),...
                    "Q",NaN(1,nFiles),...
                    "peakLocation",NaN(1,nFiles),...
                    "peakBorder",NaN(2,nFiles),...
                    "peakRTs",NaN(1,nFiles),...
                    "adductType",strings,...
                    "mass",[],...
                    "formula",strings,...
                    "signal2Noise",NaN(1,nFiles),...
                    "entropy",NaN(1,nFiles),...
                    "XIC",cell(1),...
                    "spectrumMS1",cell(1),...
                    "spectrumMS2",cell(1),...
                    "asymmetry",[]);
                EmptyStruct.spectrumMS1 = cell(1,nFiles);
                EmptyStruct.spectrumMS2 = cell(1,nFiles);
                EmptyStruct.isotopePattern = cell(1,nFiles);
                EmptyStruct.chargeState = NaN(1,nFiles);

                storedFeatures = cell(length(FeatureResults),1);

                %gather tolerances
                timeTolerance = obj.peakTimeTolerance;

                %match features and store in feature struct

                parfor iFeature = 1:length(FeatureResults)
                    currentFeatureStruct = EmptyStruct;
                    currentFeatureStruct.mz_measured = FeatureResults(iFeature).mass;
                    currentFeatureStruct.XIC = FeatureResults(iFeature).XIC;

                    nPeaks = numel(vertcat(FeatureResults(iFeature).peakLocation{:}));

                    %unpack data
                    peakData = zeros(nPeaks,10);
                    peakData(:,1) = vertcat(FeatureResults(iFeature).peakLocation{:});
                    peakData(:,2) = vertcat(FeatureResults(iFeature).peakRetentionTime{:});
                    peakData(:,3) = vertcat(FeatureResults(iFeature).peakStartLocation{:});
                    peakData(:,4) = vertcat(FeatureResults(iFeature).peakEndLocation{:});
                    peakData(:,5) = vertcat(FeatureResults(iFeature).peakHeight{:});
                    peakData(:,6) = vertcat(FeatureResults(iFeature).peakArea{:});
                    peakData(:,7) = vertcat(FeatureResults(iFeature).entropy{:});
                    peakData(:,8) = vertcat(FeatureResults(iFeature).signal2Noise{:});
                    peakData(:,9) = vertcat(FeatureResults(iFeature).fileID{:});
                    peakData(:,10) = (peakData(:,4)-peakData(:,1))./(peakData(:,1)-peakData(:,3)); %asymmetry factor

                    % preallocate current feature Storage
                    currentFeatureStruct = repmat(currentFeatureStruct,height(peakData),1);

                    %store first new entry
                    currentFile = peakData(1,9);
                    currentFeatureStruct(1).peakLocation(currentFile) = peakData(1,1);
                    currentFeatureStruct(1).peakRTs(currentFile) = peakData(1,2);
                    currentFeatureStruct(1).rt = peakData(1,2);
                    currentFeatureStruct(1).peakBorder(:,currentFile) = [peakData(1,3);peakData(1,4)];
                    currentFeatureStruct(1).Int(currentFile) = peakData(1,5);
                    currentFeatureStruct(1).Area(currentFile) = peakData(1,6);
                    currentFeatureStruct(1).entropy(currentFile) = peakData(1,7);
                    currentFeatureStruct(1).signal2Noise(currentFile) = peakData(1,8);
                    currentFeatureStruct(1).asymmetry = peakData(1,10);
                    peakData(1,:) = [];

                    %assign remaining peaks to features
                    while ~isempty(peakData)

                        currentRT = peakData(1,2);
                        currentFile = peakData(1,9);
                        currentAsymmetry = peakData(1,10);
                        id = abs(vertcat(currentFeatureStruct(:).rt) - currentRT) <= timeTolerance;
                        matchingRT = sum(id);

                        if matchingRT == 0 %no matching RT -> new Feature
                            currentFeatureStruct = storeInNewFeat(currentFeatureStruct,currentFile,peakData);

                        elseif matchingRT == 1 %single feature -> store
                            if isnan(currentFeatureStruct(id).peakLocation(currentFile)) %check if a peak is already present
                                currentFeatureStruct(id).peakLocation(currentFile) = peakData(1,1);
                                currentFeatureStruct(id).peakRTs(currentFile) = peakData(1,2);
                                currentFeatureStruct(id).peakBorder(:,currentFile) = [peakData(1,3);peakData(1,4)];
                                currentFeatureStruct(id).Int(currentFile) = peakData(1,5);
                                currentFeatureStruct(id).Area(currentFile) = peakData(1,6);
                                currentFeatureStruct(id).entropy(currentFile) = peakData(1,7);
                                currentFeatureStruct(id).signal2Noise(currentFile) = peakData(1,8);
                                %average retentionTime
                                currentFeatureStruct(id).rt = mean([currentFeatureStruct(id).rt;peakData(1,2)],'omitnan');
                            else
                                currentFeatureStruct = storeInNewFeat(currentFeatureStruct,currentFile,peakData);
                            end

                        else %multiple matching features -> store based on asymmetry factor
                            [~,idA] = min(vertcat(currentFeatureStruct(id).asymmetry) - currentAsymmetry);
                            if isnan(currentFeatureStruct(idA).peakLocation(currentFile)) %check if a peak is already present
                                currentFeatureStruct(idA).peakLocation(currentFile) = peakData(1,1);
                                currentFeatureStruct(idA).peakRTs(currentFile) = peakData(1,2);
                                currentFeatureStruct(idA).peakBorder(:,currentFile) = [peakData(1,3);peakData(1,4)];
                                currentFeatureStruct(idA).Int(currentFile) = peakData(1,5);
                                currentFeatureStruct(idA).Area(currentFile) = peakData(1,6);
                                currentFeatureStruct(idA).entropy(currentFile) = peakData(1,7);
                                currentFeatureStruct(idA).signal2Noise(currentFile) = peakData(1,8);
                                %average retentionTime
                                currentFeatureStruct(idA).rt = mean([currentFeatureStruct(idA).rt;peakData(1,2)],'omitnan');
                            else
                                currentFeatureStruct = storeInNewFeat(currentFeatureStruct,currentFile,peakData);
                            end
                        end
                        
                        %remove stored peak from list
                        peakData(1,:) = [];
                    end

                    %remove empty struct
                    id = isnan([currentFeatureStruct(:).rt])';
                    currentFeatureStruct(id) = [];
                    % get final feature rt
                    for jStruct = 1:numel(currentFeatureStruct)
                        currentFeatureStruct(jStruct).rt = mean(currentFeatureStruct(jStruct).peakRTs,"all","omitmissing");
                    end
                    %store currentFeatureStruct
                    storedFeatures{iFeature,1} = currentFeatureStruct;
                end

                %unzip features
                storedFeatures = vertcat(storedFeatures{:});

                if isISIntegration == false
                    %gather original scans
                    storedFeatures = obj.findOriginalMassScans(storedFeatures);
                    % trim XIC to relevant parts
                    storedFeatures = obj.trimExtractedIonChromatograms(storedFeatures);
                    % extract isotope distribution
                    storedFeatures = obj.gatherIsotopeDistributions(storedFeatures);
                    % confirm same feature
                    storedFeatures = obj.confirmSameFeatureByIsotopeDistribution(storedFeatures);
                    % Occurrence filter
                    [storedFeatures,obj.occurenceFiltered] = obj.occurrenceFilterFeatures(storedFeatures);
                    % calculate Peak width MZ dimension
                    storedFeatures = obj.calculatePeakWidth_MZ(storedFeatures);
                    % calculate Peak width Time dimension
                    storedFeatures = obj.calculatePeakWidth_Time(storedFeatures);
                    % build feature isotope pattern
                    storedFeatures = obj.averageIsotopePattern(storedFeatures);
                    storedFeatures = obj.correctMassByChargeState(storedFeatures);
                    % gather MS2 scans
                    storedFeatures = obj.gatherMS2Spectra(storedFeatures);
                end

                FeatureResults = storedFeatures;
            else
                FeatureResults = IntegrationResults;
            end
        end

        %% Hard Ionisation Functions
        function [Output,obj] = extractFeatures_HardIonisation(obj,progressBar)


            if obj.useMassAlign == true
                progressBar.Message = "Aligning MS Scans";
                obj = obj.alignMasses("batch");
                progressBar.Value = 0.4;
            end

            % ROI Search
            progressBar.Message = "Searching for ROIs";
            obj = obj.findRegionOfInterest("batch");
            progressBar.Value = 0.5;

            % Baseline Correction
            if obj.useBaselineCorrection == true
                progressBar.Message = "Correcting Baseline";
                obj = obj.correctBaseline("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end
            % Smoothing
            if obj.useSmoothing == true
                progressBar.Message = "Smoothing Peaks";
                obj = obj.smoothPeaks("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Peak Align
            if obj.usePeakAlign == true && nData > 1
                progressBar.Message = "Aligning Peaks";
                obj = obj.alignPeaks("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Blank Correction
            if any(contains(obj.fileType,"blank")) % Separate Blank data from Sample data and subtract
                tempBLK = obj.TempDataFileObj.ROICells(end,1);
                obj.TempDataFileObj.ROIMatBLK = sparse(tempBLK{:});
                obj.TempDataFileObj.ROICells(end) = [];
                obj.TempDataFileObj.TimeCells(end) = [];
                progressBar.Message = "Subtracting Blank";
                peakCells = obj.TempDataFileObj.ROICells;
                blankMat = obj.TempDataFileObj.ROIMatBLK;
                parfor iFile = 1:size(peakCells,1)
                    peakCells{iFile,1} = peakCells{iFile,1}-blankMat;
                    % set possible negative values to 0
                    peakCells{iFile,1} = max(peakCells{iFile,1},0);
                end
                obj.TempDataFileObj.ROICells = peakCells;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % pad arrays with Maximum peak width*3 Scans to eliminate
            % integration interference between matrices
            obj = obj.finalizeROI;

            clearvars -except obj progressBar
            %% Integration Stage
            % Find and Integrate IS separate
            if obj.useInternalStandard == true
                progressBar.Message = "Searching for Internal Standard";
                obj = obj.identifyInternalStandard;
                if ~isempty(obj.internalStandardData)
                    obj = obj.internalStandardNormalization;
                end
                progressBar.Value = progressBar.Value + 0.05;
            end
            
            % mass correction
            if obj.useInternalStandard == true && ~isempty(obj.interalStandardIntensity)
                progressBar.Message = "Performing IS mass correction";
                obj = obj.massCorrectionByInternalStandard;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Integrate all Peaks
            progressBar.Message = "Picking Peaks";
            IntegrationData = obj.findPeaks_Hard;
            progressBar.Value = 0.9;

            progressBar.Message = "Processing found Features";
            IntegrationData = obj.gatherEISpectra(IntegrationData);
            IntegrationData = obj.assignFileID(IntegrationData);
            IntegrationData = obj.fileSortPeaks(IntegrationData);
            IntegrationData = obj.mergeDuplicatePeaksWithinFile(IntegrationData);

            Output = obj.initializeOutputStruct(IntegrationData);
            % Build Storage Arrays
            FeatureData = obj.buildFeatureArray_Hard(IntegrationData);

            % confirm same feature by MS2 comparison
            FeatureData = obj.confirmSameFeatureByMS2(FeatureData);

            % Occurrence filter
            [FeatureData,obj.occurenceFiltered] = obj.occurrenceFilterFeatures(FeatureData);

            %build average EI (MS2) spectrum
            FeatureData = obj.finalizeEISpectra(FeatureData);
            progressBar.Value = 0.95;

            %apply scaling
            if obj.useScaling == true
                progressBar.Message = "Apply scaling";
                FeatureData = obj.groupAndSampleScaling(FeatureData);
            end

            %fill remaining fields
            Output = obj.finalizeOutputStruct(Output,FeatureData);
            progressBar.Message = "Group processing successful";
            progressBar.Value = 1;

            obj.Output = Output;
            delete(progressBar)
        end

        function IntegrationResults = findPeaks_Hard(obj)
            %output preallocation
            IntegrationResults = struct( ...
                "mass",[], ...
                "peakLocation",[], ...
                "peakRetentionTime",[], ...
                "peakStartLocation",[], ...
                "peakEndLocation",[], ...
                "peakHeight",[], ...
                "peakArea",[], ...
                "entropy",[], ...
                "signal2Noise",[], ...
                "minWidthFiltered",[], ...
                "maxWidthFiltered",[], ...
                "entropyFiltered",[], ...
                "signal2NoiseFiltered",[], ...
                "spectrumMS2",[], ...
                "XIC",[], ...
                "fileID",[]);

            %prepare TIC Data
            tics = sum(obj.TempDataFileObj.ROIMat,2);
            times = obj.TempDataFileObj.timeVec;
            IntegrationResults.XIC = [times,tics];
            %gather parameters
            currentTIC = full(tics);
            currentTime = full(times);
            smoothedTIC = smoothdata(currentTIC,"gaussian",4,"omitnan");
            minPeakWidthScans = floor(obj.peakMinWidth/obj.scanFrequencySecond);
            %calculate noise level (10th percentile of non zero values)
            noise = prctile(currentTIC(currentTIC > 0),10);
            [~,peakLoc,peakWidth] = findpeaks(currentTIC,"WidthReference","halfheight");
            %calculate initial borders and bring in correct form
            lowerBorders = max(floor(peakLoc-peakWidth/2),1); % limit lower peak border to scan 1
            upperBorders = min(ceil(peakLoc+peakWidth/2),numel(currentTIC)); % limit upper peak border to last scan
            peakData = [peakLoc,lowerBorders,upperBorders];
            peakData = correctPeakData(peakData,currentTIC,smoothedTIC,minPeakWidthScans);
            IntegrationResults.peakLocation = peakData(:,1);
            IntegrationResults.peakStartLocation = peakData(:,2);
            IntegrationResults.peakEndLocation = peakData(:,3);
            IntegrationResults.peakHeight = peakData(:,4);
            IntegrationResults = obj.filterPeaksFromIntegration(IntegrationResults,noise);
            % entropy calculation
            IntegrationResults = obj.calculatePeakEntropy(IntegrationResults);
            IntegrationResults = obj.integratePeaks(IntegrationResults,currentTime);
        end

        function FeatureData = buildFeatureArray_Hard(obj,IntegrationResults,varargin)
            nFiles = numel(obj.fileName);

            xicData = IntegrationResults(1).XIC;

            %remove unnecessary fields from input struct
            IntegrationResults = rmfield(IntegrationResults,["minWidthFiltered","maxWidthFiltered","entropyFiltered","signal2NoiseFiltered","XIC"]);

            featureStruct = struct(...
                "featID",strings,...
                "mass_measured",[],...
                "retentionTime",[],...
                "adductType",strings,...
                "mass_corrected",[],...
                "formula",strings,...
                "peakHeights",zeros(0,nFiles),...
                "peakAreas",zeros(0,nFiles),...
                "peakLocations",zeros(0,nFiles),...
                "peakBorders",zeros(0,nFiles),...
                "retentionTimes",zeros(0,nFiles),...
                "signal2Noise",zeros(0,nFiles),...
                "entropy",zeros(0,nFiles),...
                "XIC",cell(1),...
                "spectrumMS1",cell(1),...
                "spectrumMS2",cell(1));

            %gather tolerances
            TIMETOLERANCE = obj.peakTimeTolerance;
            MASSTOLERANCE = obj.betweenFileMassTolerance;
            MASSUNIT = obj.betweenFileMassUnit;

            if isscalar(varargin)
                isISIntegration = true;
            else
                isISIntegration = false;
            end

            %match features and store in feature struct

            uniqueFeatures = [vertcat(IntegrationResults(:).mass),vertcat(IntegrationResults(:).peakRetentionTime)];

            % over preallocat feature Storage
            featureStruct = repmat(featureStruct,height(uniqueFeatures),1);
            counter = 0;
            while ~isempty(uniqueFeatures)
                counter = counter + 1;
                currentFeature = uniqueFeatures(1,:);
                % preallocate temp storage
                emptyArray = NaN(1,nFiles);
                areas = emptyArray;
                heights = emptyArray;
                retentionTimes = emptyArray;
                peakLocation = emptyArray;
                peakBorders = [emptyArray;emptyArray];
                signal2Noise = emptyArray;
                entropy = emptyArray;
                spectrum = cell(1,nFiles);
                xic = cell(1,nFiles);

                %compare feature between files
                for iFile = 1:nFiles
                    switch MASSUNIT
                        case "Da"
                            idm = abs(IntegrationResults(iFile).mass-currentFeature(1,1)) <= MASSTOLERANCE;
                        case "ppm"
                            idm = abs(IntegrationResults(iFile).mass-currentFeature(1,1))./currentFeature(1,1)*10^6 <= MASSTOLERANCE;
                    end
                    idRT = abs(IntegrationResults(iFile).peakRetentionTime - currentFeature(1,2)) <= TIMETOLERANCE;
                    idx = idm & idRT;
                    % handle matching peaks
                    if sum(idx) == 0 %no matching peaks
                        continue

                    elseif sum(idx) == 1 && ~all(isnan(peakLocation)) %compare spectra if already peak stored in feature
                        %build current average spectrum
                        currentSpectrum = alignSpectra(spectrum,"average","low","true");
                        %gather possible spectra
                        possibleSpectra = IntegrationResults(iFile).spectrumMS2(idx);
                        possibleSpectra = alignSpectra(possibleSpectra,"normal","low","true");
                        %calculate scores
                        [CompoundScores,~] = scoresBetweenSets(currentSpectrum,possibleSpectra);
                        if all(CompoundScores(:,1) < 700) %features don´t match
                            continue
                        end

                    elseif sum(idx) > 1 %use feat with higher Similarity score
                        %build current average spectrum
                        currentSpectrum = alignSpectra(spectrum,"average","low","true");
                        %gather possible spectra
                        possibleSpectra = IntegrationResults(iFile).spectrumMS2(idx);
                        possibleSpectra = alignSpectra(possibleSpectra,"normal","low","true");
                        %calculate scores
                        [CompoundScores,~] = scoresBetweenSets(currentSpectrum,possibleSpectra);
                        if all(CompoundScores(:,1) < 700) %no matching feature
                            continue
                        end
                        %get id of feature with higher score
                        idmax = CompoundScores(CompoundScores(:,1)==max(CompoundScores(:,1)),3);
                        location = find(idx);
                        idx = false(size(idx));
                        idx(location(idmax)) = true;
                    end
                    %store found peak information
                    areas(1,iFile) = IntegrationResults(iFile).peakArea(idx);
                    heights(1,iFile) = IntegrationResults(iFile).peakHeight(idx);
                    retentionTimes(1,iFile) = IntegrationResults(iFile).peakRetentionTime(idx);
                    peakLocation(1,iFile) = IntegrationResults(iFile).peakLocation(idx);
                    peakBorders(1,iFile) = IntegrationResults(iFile).peakStartLocation(idx);
                    peakBorders(2,iFile) = IntegrationResults(iFile).peakEndLocation(idx);
                    signal2Noise(1,iFile) = IntegrationResults(iFile).signal2Noise(idx);
                    entropy(1,iFile) = IntegrationResults(iFile).entropy(idx);
                    spectrum(1,iFile) = IntegrationResults(iFile).spectrumMS2(idx);
                    xic{1,iFile} = full(xicData(peakBorders(1,iFile):peakBorders(2,iFile),:));

                    %delete peaks from input struct
                    IntegrationResults(iFile).mass(idx) = [];
                    IntegrationResults(iFile).peakArea(idx) = [];
                    IntegrationResults(iFile).peakHeight(idx) = [];
                    IntegrationResults(iFile).peakRetentionTime(idx) = [];
                    IntegrationResults(iFile).peakLocation(idx) = [];
                    IntegrationResults(iFile).peakStartLocation(idx) = [];
                    IntegrationResults(iFile).peakEndLocation(idx) = [];
                    IntegrationResults(iFile).signal2Noise(idx) = [];
                    IntegrationResults(iFile).entropy(idx) = [];
                    IntegrationResults(iFile).spectrumMS2(idx) = [];


                end
                %store matching features
                featureStruct(counter).mass_measured = currentFeature(1,1);
                featureStruct(counter).retentionTime = currentFeature(1,2);
                featureStruct(counter).peakHeights = heights;
                featureStruct(counter).peakAreas = areas;
                featureStruct(counter).peakLocations = peakLocation;
                featureStruct(counter).peakBorders = peakBorders;
                featureStruct(counter).retentionTimes = retentionTimes;
                featureStruct(counter).signal2Noise = signal2Noise;
                featureStruct(counter).entropy = entropy;
                featureStruct(counter).spectrumMS2 = spectrum;
                featureStruct(counter).XIC = xic;

                %update remaining features
                uniqueFeatures(1,:) = [];
            end

            if isISIntegration == false

                %gather original scans
                featureStruct = obj.findOriginalMassScans(featureStruct);
            end

            FeatureData = featureStruct;
        end
        %% Furter Processing
        function outputStruct = gatherMS2Spectra(obj,FeatureStruct)
            %check if MSn data exists, skip if not
            if isscalar(obj.RawDataFileObj.centroidDataMS2)
                outputStruct = FeatureStruct;
                return
            end

            %%%%%
            % test tolerances
            MASSTOLERANCE = 0.05;
            MASSUNIT = "Da";
            TIMETOLERANCE = 10;
            %%%%%%

            timeArray = obj.RawDataFileObj.timeDataMS2;
            timeArray = vertcat(timeArray{:});
            scanArray = obj.RawDataFileObj.centroidDataMS2;
            scanArray = vertcat(scanArray{:});
            precursor = obj.RawDataFileObj.molecularPrecursorMass;
            precursor = vertcat(precursor{:});

            parfor iFeature = 1:length(FeatureStruct)
                idMass = [];
                switch MASSUNIT
                    case "Da"
                        idMass = abs(precursor-FeatureStruct(iFeature).mass_measured) <= MASSTOLERANCE;
                    case "ppm"
                        idMass = abs(precursor-FeatureStruct(iFeature).mass_measured)./FeatureStruct(iFeature).mass_measured*10^6 <= MASSTOLERANCE;
                end
                idTime = abs(timeArray-FeatureStruct(iFeature).retentionTime) <= TIMETOLERANCE;
                id = idTime & idMass;
                foundScan = scanArray(id);
                %remove possible empty scans
                foundScan(cellfun(@isempty, foundScan)) = [];
                if numel(foundScan) >= 1
                    %remove masses > precursor mass
                    for jScan = 1:height(foundScan)
                        highMassId = foundScan{jScan,1}(:,1) > FeatureStruct(iFeature).mass_measured + 1;
                        foundScan{jScan,1}(highMassId,:) = [];
                    end
                    foundScan = alignSpectra(foundScan,"average","low","true",0.05);
                else % no found scan
                    foundScan = [];
                end
                FeatureStruct(iFeature).spectrumMS2 = foundScan;
            end
            outputStruct = FeatureStruct;
        end

        function [IntResults,obj] = gatherEISpectra(obj,IntResults)
            intensities = obj.TempDataFileObj.ROIMat;
            masses =  obj.TempDataFileObj.ROImzVec;
            foundSpectrumArray = cell(size(IntResults.peakLocation));
            molecularMass = zeros(size(foundSpectrumArray));
            peakWidths = [IntResults.peakStartLocation,IntResults.peakEndLocation];
            EIlosses = load("MassListData.mat","EICommonLoss");
            EIlosses = EIlosses.EICommonLoss;
            parfor iMass = 1:height(molecularMass)
                peakBorders = peakWidths(iMass,:);
                spectrum = intensities(peakBorders(1):peakBorders(2),:);
                %normalize Spectras
                spectrum = spectrum./max(spectrum,[],"all");
                %remove columns with more than 50% empty
                id = (sum(spectrum~=0)/height(spectrum))<0.5
                spectrum(:,id) = []
                %mean spectra
                spectrum = mean(spectrum);
                spectrum = [masses(~id);full(spectrum)]';
                %remove rows with intensity < 0.01
                spectrum(spectrum(:,2)<0.01,:) = [];
                foundSpectrumArray{iMass,1} = spectrum;

                %identify molecular mass
                spectrumMass = flip(spectrum(:,1));
                hasMolecularMass = false;
                counter = 0;
                while hasMolecularMass == false & counter < numel(spectrumMass)
                    counter = counter + 1;
                    possibleFragment = spectrumMass(counter)-EIlosses;
                    hasMolecularMass = any(min(abs(spectrumMass-possibleFragment'))<0.1);
                end
                if hasMolecularMass
                    molecularMass(iMass,1) = spectrumMass(counter,1);
                else
                    molecularMass(iMass,1) = spectrum(end,1);
                end
            end
            IntResults.spectrumMS2 = foundSpectrumArray;
            IntResults.mass = round(molecularMass,1);
        end

        function OutArray = mergeDuplicatePeaksWithinFile(obj,InArray)

            MASSTOLERANCE = obj.withinFileMassTolerance;
            MASSUNIT = obj.withinFileMassUnit;
            TIMETOLERANCE = obj.peakTimeTolerance;

            nFeature = size(InArray.mass);

            OutArray = struct( ...
                "mass",cell(nFeature), ...
                "peakLocation",cell(nFeature), ...
                "peakRetentionTime",cell(nFeature), ...
                "peakStartLocation",cell(nFeature), ...
                "peakEndLocation",cell(nFeature), ...
                "peakHeight",cell(nFeature), ...
                "peakArea",cell(nFeature), ...
                "entropy",cell(nFeature), ...
                "signal2Noise",cell(nFeature), ...
                "minWidthFiltered",InArray.minWidthFiltered, ...
                "maxWidthFiltered",InArray.maxWidthFiltered, ...
                "entropyFiltered",InArray.entropyFiltered, ...
                "signal2NoiseFiltered",InArray.signal2NoiseFiltered, ...
                "spectrumMS2",cell(nFeature), ...
                "XIC",InArray.XIC, ...
                "fileID",cell(nFeature));

            fnames = fieldnames(OutArray);
            fnames(any(fnames == ["XIC","minWidthFiltered","maxWidthFiltered","entropyFiltered","signal2NoiseFiltered"],2)) = []; %remove names from list to skip field in assignmelt loop

            xic = full(InArray.XIC(:,1));

            for iFeature = 1:numel(InArray.mass) %sample loop
                featureID = [InArray.mass{1,iFeature},InArray.peakRetentionTime{1,iFeature}];
                counter = 0;
                while ~isempty(featureID)
                    counter = counter + 1;
                    % find all features that match the current signiture
                    %mass tolerance
                    switch MASSUNIT
                        case "Da"
                            idMass = abs(featureID(:,1)-featureID(1,1)) <= MASSTOLERANCE;
                        case "ppm"
                            idMass = abs(featureID(:,1)-featureID(1,1))./featureID(1,1)*10^6 <= MASSTOLERANCE;
                    end
                    %time tolerance
                    idTime = abs(featureID(:,2)-featureID(1,2)) <= TIMETOLERANCE;
                    id = idMass & idTime;

                    %% multiple peaks found
                    % remove peaks with peak height < 3x baseline
                    if sum(id) > 1
                        bordersStart = InArray.peakStartLocation{1,iFeature};
                        bordersEnd = InArray.peakEndLocation{1,iFeature};
                        %remove less prominent peak
                        idx = zeros(size(id));
                        for jPeak = 1:numel(id)
                            if id(jPeak) == false
                                continue
                            else
                                tempTIC = xic(bordersStart(jPeak):bordersEnd(jPeak),:);
                                idx(jPeak) =  mean([tempTIC(1);tempTIC(end)]) / max(tempTIC);
                            end
                        end
                        %remove feat and return while loop
                        id = id & idx(idx==max(idx));
                        for jFieldName = 1:numel(fnames) % loop over each field name
                            InArray.(fnames{jFieldName}){1,iFeature}(id) =  [];
                        end

                        %% only one peak remaining
                        % store in output and remove from input
                    elseif sum(id) == 1
                        for jFieldName = 1:numel(fnames) % loop over each field name
                            OutArray(iFeature).(fnames{jFieldName})= vertcat(OutArray(iFeature).(fnames{jFieldName}),InArray.(fnames{jFieldName}){1,iFeature}(id));
                            InArray.(fnames{jFieldName}){1,iFeature}(id) =  [];
                        end
                    else %no matching peak

                    end
                    % remove entries from current list
                    featureID(id,:) = [];
                end
            end
        end

        function FeatureStruct = correctFeatureMass(FeatureStruct) %%%% WIP %%%%
            %% WIP
            parfor iFeature = 1:height(FeatureStruct)
                FeatureStruct(iFeature).mass_corrected = FeatureStruct(iFeature).mass_measured * FeatureStruct(iFeature).chargeState;
            end
        end

        function obj = blankCorrection(obj)
            %% blankCorrection performs automatic blank subtraction for full sequences
            % Takes reapeating blank files in the TempDataFile, calculates the
            % average and replaces each block with the average blank.

            %gather data
            dataType = obj.fileType;
            roiCells = obj.TempDataFileObj.ROICells;
            timeCells = obj.TempDataFileObj.TimeCells;
            %% subtract system Blanks (all files)
            isSystemBlank = strcmp(dataType,"suitability blank");
            nBlank = sum(isSystemBlank);
            blankCells = roiCells(isSystemBlank);
            roiCells = roiCells(~isSystemBlank);
            blankMat = zeros(size(blankCells{1}));
            if nBlank > 1 % average blank files
                for iBlank = 1:numel(blankCells)
                    blankMat = blankMat + blankCells{iBlank,1};
                end
                blankMat = blankMat / nBlank;

            else  % extract cell
                blankMat = blankCells{:};
            end

            for iFile = 1:height(roiCells)
                roiCells{iFile,1} = max(roiCells{iFile,1} - blankMat,0);
            end

            % save 
            obj.TempDataFileObj.ROIMatSystemBLK = sparse(blankMat);
            % remove processed blank data
            timeCells(isSystemBlank) = []; 
            obj.nScan(isSystemBlank) = [];

            % flag files as non sample
            obj.isSampleFile(isSystemBlank) = false;
            %% subtract process blanks from segments between blanks

            % Find indices of "process blank"
            blankIdx = find(strcmp(dataType, "process blank"));
            maxIndex = numel(roiCells);

            for iBlank = 1:numel(blankIdx)
                idxBlank = blankIdx(iBlank);

                % start just after this blank
                startIdx = idxBlank + 1;

                % end before the next blank, or to the end if this is the last one
                if iBlank < numel(blankIdx)
                    endIdx = blankIdx(iBlank+1) - 1;
                else
                    endIdx = maxIndex;
                end

                % Value of this process blank
                blankMat = roiCells{idxBlank};

                % If nothing after this blank, skip
                if startIdx > endIdx
                    continue
                else

                % Subtract from all following cells in this block
                for jFile = startIdx:endIdx
                    roiCells{jFile} = max(roiCells{jFile} - blankMat,0);
                end
                end
            end
            obj.TempDataFileObj.ROIMatProcessBLK = roiCells(blankIdx);
            
            %remove processed blank data from obj
            timeCells(blankIdx) = [];
            roiCells(blankIdx) = [];
            obj.nScan(blankIdx) = [];
            % flag files as non sample
            obj.isSampleFile(blankIdx) = false;

            % store corrected data on disk
            obj.TempDataFileObj.TimeCells = timeCells;
            obj.TempDataFileObj.ROICells = roiCells;
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

        function obj = internalStandardNormalization(obj) %%%% MUST BE UPDATED %%%%
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

        function obj = identifyInternalStandard(obj) %%%% MUST BE UPDATED %%%%
            % identify IS Vectors
            obj.ISMass = obj.ISDat(:,1)';
            mzVec = obj.TempDataFileObj.ROImzVec;
            %check tolerance and store difference
            switch obj.mzTolUnit
                case "Da"
                    ISid = abs(mzVec-obj.ISMass') <= obj.mzTol;
                    obj.ISMassFound = sort(mzVec(any(ISid,1)),2,"ascend");
                    foundMassID = any(ISid,2);
                case "ppm"
                    ISid = abs(mzVec-obj.ISMass')./mzVec*10^6 <= obj.mzTol;
                    obj.ISMassFound = sort(mzVec(any(ISid,1)),2,"ascend");
                    foundMassID = any(ISid,2);
            end
            %remove IS outside tolerance and throw warning
            if all(~foundMassID)
                % if no IS mass found, throw warning and exit
                header="Skipping ISTD Normalization";
                message = ["Reason: no Internal standard mass was found in this group","Check the specified m/z or increase the mass tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
                return
            elseif any(~foundMassID)
                % if some IS mass is not found, throw warning and continue
                header="Skipping ISTD No. ";
                for i=1:size(obj.ISMass,2)
                    if foundMassID(i)==0
                        header = header + i + " ";
                    end
                end
                message = ["Reason: Internal standard mass was not found in this group","Check the specified m/z or increase the mass tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
            end
            ISid=any(ISid);
            obj.ISMass(~foundMassID)=[];
            % extract relevant columns and perform Peak Picking and
            % Integration
            ISIntegrationData = obj.findPeaks(ISid);
            % remove possible empty columns
            id = cellfun(@isempty,ISIntegrationData(1,:));
            ISIntegrationData(:,id) = [];
            obj.ISMassFound(:,id) = [];
            if isempty(ISIntegrationData{1,1})
                header="Skipping ISTD Normalization";
                message = ["Reason: no Internal standard was found in this group","Check the specified RT or increase the time tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
                return
            end
            % assign Peaks to Sample
            ISIntegrationData = obj.AssignRT2SampleFile(ISIntegrationData);
            ISIntegrationData(5:7,:) = [];
            ISData = obj.buildFeatureArray_Soft(ISIntegrationData,obj.ISMassFound);

            % check if RT Range is Correct and Remove Feature outside range
            nIS = height(obj.ISDat);
            counter = 1;
            id = [];
            while counter <= nIS
                ISmz = obj.ISDat(counter,1);
                IStime = obj.ISDat(counter,2);
                timeTol = obj.ISDat(counter,3);
                switch obj.mzTolUnit
                    case "Da"
                        idmz = abs(ISData.FeatIdentifiers(:,1)-ISmz) >= obj.mzTol;
                    case "ppm"
                        idmz = abs(ISData.FeatIdentifiers(:,1)-ISmz)./ISmz*10^6 >= obj.mzTol;
                end
                idrt = abs(ISData.FeatIdentifiers(:,2)-IStime)>=timeTol;
                id = [id,any([idmz,idrt],2)];
                counter = counter+1;
            end
            id = all(id,2);
            ISData.FeatIdentifiers(id,:) = [];
            ISData.IntensityStorage(id,:) = [];
            ISData.RetentionTimeStorage(id,:) = [];
            ma = [];
            % filter possible multiple Features for one mass
            if ~isempty(ISData.FeatIdentifiers)
                meanInt = mean(ISData.IntensityStorage,2);
                counter = 1;
                while counter <= nIS
                    ISmz = obj.ISDat(counter,1);
                    switch obj.mzTolUnit
                        case "Da"
                            idmz = abs(ISData.FeatIdentifiers(:,1)-ISmz) <= obj.mzTol;
                        case "ppm"
                            idmz = abs(ISData.FeatIdentifiers(:,1)-ISmz)./ISData.FeatIdentifiers(:,1)*10^6 <= obj.mzTol;
                    end
                    ma = [ma,max(meanInt(idmz))];
                    counter = counter+1;
                end
                id = ismember(meanInt,ma);
                ISData.FeatIdentifiers = ISData.FeatIdentifiers(id,:);
                ISData.IntensityStorage = ISData.IntensityStorage(id,:);
                ISData.RetentionTimeStorage = ISData.RetentionTimeStorage(id,:);
            end

            %store in obj
            check = ismember(obj.ISMassFound,ISData.FeatIdentifiers(:,1));
            obj.ISMassFound(~check) = [];
            obj.ISRT = ISData.FeatIdentifiers(:,2)';
            obj.ISValue = ISData.IntensityStorage';

            %check wich IS remains
            foundMassID = false(size(obj.ISMass));
            for n = 1:width(obj.ISMass)
                ISmz = obj.ISMass(1,n);
                switch obj.mzTolUnit
                    case "Da"
                        foundMassID(n) = any(abs(obj.ISMassFound-ISmz) <= obj.mzTol);
                    case "ppm"
                        foundMassID(n) = any(abs(obj.ISMassFound-ISmz)./ISmz*10^6 <= obj.mzTol);
                end
            end
            obj.ISMass = sort(obj.ISMass(foundMassID),2,'ascend');
            obj.ISdelta = obj.ISMass-obj.ISMassFound;
            if all(~foundMassID)
                % if no IS mass found, throw warning and exit
                header="Skipping ISTD Normalization";
                message = ["Reason: Every internal standard is missing peaks in one or more samples.","Check the specified retention time or increase the time tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
                return
            elseif any(~foundMassID)
                % if some IS mass is not found, throw warning and continue
                header="Skipping ISTD No. ";
                for i=1:size(foundMassID,2)
                    if foundMassID(i) == 0
                        header = header + i + ", ";
                    end
                end
                message = ["Reason: Internal standard is missing peaks in one or more samples","Check the specified retention time or increase the time tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
            end
        end

        function obj = massCorrectionByInternalStandard(obj) %%%% MUST BE UPDATED %%%%
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

        function [valuesFiltered,obj] = filterAdducts(obj,IntegrationResults) %%%% MUST BE UPDATED %%%%
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
            tempPeakData = obj.TempDataFileObj.ROICells;
            tempTimeData = obj.TempDataFileObj.TimeCells;

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
                outROI{1,1} = smoothdata(msRoi_end,"gaussian",3); %apply slight smoothing, to remove gaps within peaks
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
            QUANTIL = obj.baselineQuantil;

            parfor iFile = 1:size(msRoi,1)
                oldSize = height(msRoi{iFile,1});
                %depad Array
                msRoiTemp = msRoi{iFile,1};
                [msRoiTemp,timeTemp] = depadArrays(msRoiTemp,time{iFile,1});
                msRoiTemp = msbackadj(timeTemp,msRoiTemp,'WindowSize',WINDOW_SIZE,'StepSize',STEP_SIZE,'RegressionMethod',REGRESSION,'EstimationMethod',ESTIMATION,'SmoothMethod','none','QuantileValue',QUANTIL,'PreserveHeights',true);
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

        function OutputStruct = groupAndSampleScaling(obj,FeatureStruct)
            %% applies group and sample scaling to final output struct

            GROUP_SCALING_FACTORS = obj.groupScale;
            SAMPLE_SCALING_FACTORS = obj.sampleScale';

            if InputStruct.dataSize > 0
                parfor iFeature = 1:length(FeatureStruct)
                    %GroupScale
                    FeatureStruct(iFeature).peakHeights = FeatureStruct(iFeature).peakHeights/GROUP_SCALING_FACTORS;
                    FeatureStruct(iFeature).peakAreas = FeatureStruct(iFeature).peakAreas/GROUP_SCALING_FACTORS;
                    %SampleScale
                    FeatureStruct(iFeature).peakHeights = FeatureStruct(iFeature).peakHeights./SAMPLE_SCALING_FACTORS;
                    FeatureStruct(iFeature).peakAreas = FeatureStruct(iFeature).peakAreas./SAMPLE_SCALING_FACTORS;
                end
                OutputStruct = FeatureStruct;
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

                %remove peaks with height = 1 (no peak location was found)
                idToRemoveNew = IntResults(iFeature).peakHeight == 1;
                idxToRemove = idxToRemove | idToRemoveNew;

                % %remove peaks with bad Peak asymmetry
                % symmetry = (IntResults(iFeature).peakEndLocation - IntResults(iFeature).peakLocation)./(IntResults(iFeature).peakLocation - IntResults(iFeature).peakStartLocation);
                % idToRemoveNew = symmetry<0.3 | symmetry>3;
                % idxToRemove = idxToRemove | idToRemoveNew;

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

        function OutputStruct = findOriginalMassScans(obj,FeatureStruct)
            %% gathers MS1 spectra for each feature
            OutputStruct = FeatureStruct;
            % load data
            allDataStruct = obj.RawDataFileObj;
            % preallocation
            nFile = numel(allDataStruct);
            allScans = cell(nFile,1);
            % append all scans with spacers in between, to match processing
            % indices
            for iScan = 1:nFile
                temp = {allDataStruct(iScan).spectraMS1.processedScan}';
                temp(obj.nScanPadded(iScan),1) = {[]};
                allScans{iScan,1} = temp;
            end
            allScans = vertcat(allScans{:});
            
            %% find original MS1 spectra
            parfor iFeature = 1:numel(FeatureStruct)
                spectra = cell(1,nFile);
                location = FeatureStruct(iFeature).peakBorder;

                for jFile = 1:nFile
                    %check if borders contain NaN then skip iteration
                    if any(isnan(location(:,jFile)))
                        continue
                    end
                    %select spectra in peak range
                    scans = allScans(location(1,jFile):location(2,jFile));
                    %remove possible empty scans
                    scans(cellfun(@isempty, scans)) = [];
                    if ~isempty(scans)
                        spectra{1,jFile} = scans;
                    else
                        spectra{1,jFile} = {[]};
                    end
                end
                OutputStruct(iFeature).spectrumMS1 = spectra;
            end
        end

        function [Output,sumFiltered] = occurrenceFilterFeatures(obj,FeatureStruct)
            %% remove features with less peaks than required minimum from Output struct
            nFile = numel(obj.dataFile);
            minDataPoints = ceil(nFile*obj.minOccurence);
            numElements = zeros(length(FeatureStruct),1);

            for iFeature = 1:length(FeatureStruct)
                numElements(iFeature) = nnz(~isnan(FeatureStruct(iFeature).peakHeights));
            end
            idToRemove = numElements < minDataPoints;

            %sum number of removed peaks
            sumFiltered = sum(numElements(idToRemove),"all");
            FeatureStruct(idToRemove) = [];

            Output = FeatureStruct;
        end

        function Output = initializeOutputStruct(obj,IntegrationData)
            %preallocate Output struct
            Output = struct(...
                "feature",[],...
                "minWidthFiltered",[],...
                "maxWidthFiltered",[],...
                "entropyFiltered",[],...
                "signal2NoiseFiltered",[],...
                "occurrenceFiltered",[],...
                "groupName",string,...
                "fileNames",string,...
                "dataSize",[],...
                "separationType",string);

            %store group infos
            Output.minWidthFiltered = sum(vertcat(IntegrationData(:).minWidthFiltered));
            Output.maxWidthFiltered = sum(vertcat(IntegrationData(:).maxWidthFiltered));
            Output.entropyFiltered = sum(vertcat(IntegrationData(:).entropyFiltered));
            Output.signal2NoiseFiltered = sum(vertcat(IntegrationData(:).signal2NoiseFiltered));

            Output.fileNames = obj.fileName;
            Output.groupName = obj.groupName;
            Output.separationType = obj.separationType;

        end

        function Output = finalizeOutputStruct(obj,Output,featureStruct)
            %% finalizes feature Output Struct
            % fills remaining fields: dataSize, featID
            % rounds retentionTime
            % not yet implemented: formula, adduct type, corrected mass

            Output.dataSize = length(featureStruct);
            Output.occurrenceFiltered = obj.occurenceFiltered;

            %build average RetentionTime and featureID
            for iFeat = 1:Output.dataSize
                %build retention time
                featureStruct(iFeat).retentionTime = mean(featureStruct(iFeat).retentionTimes,"all","omitmissing");
                %build feature ID string
                featureStruct(iFeat).featID = featureStruct(iFeat).mass_measured + "Da@" + featureStruct(iFeat).retentionTime + "s_" + obj.groupName;
            end

            % separation Type specific tasks

            switch Output.separationType
                case "GC"
                    for iFeat = 1:Output.dataSize
                        % calculate formula
                        % GC -> database Search

                        % add Adduct type
                        featureStruct(iFeat).adductType = "M+";
                        % calculate corrected mass based on adduct Type or Database
                    end

                otherwise %ESI
                    
                    % determine Adduct type

                    % calculate corrected mass based on adduct Type or Database

                    % calculate formula
                    featureStruct = obj.calculateFormulaFromMass(featureStruct);

                    
                    % ESI check ms1 spectrum

                    
            end
            %store final struct
            Output.feature = featureStruct;

        end

        function FeatureStruct = gatherIsotopeDistributions(obj,FeatureStruct)

            TOLERANCE = obj.withinFileMassTolerance;
            TOLERANCEUNIT = obj.withinFileMassUnit;

            parfor iFeature = 1:height(FeatureStruct)
                targetMass = FeatureStruct(iFeature).mass_measured;

                if strcmp(TOLERANCEUNIT,"ppm")
                    adjustedTolerance = TOLERANCE * targetMass * 10^-6;
                else
                    adjustedTolerance = TOLERANCE;
                end

                spectrum = FeatureStruct(iFeature).spectrumMS1;
                isotopePattern = cell(1,numel(spectrum));
                chargeState = zeros(1,numel(spectrum));
                for jFile = 1:numel(spectrum)
                    if ~isempty(spectrum{1,jFile})
                        [isotopePattern{1,jFile},chargeState(1,jFile)] = extractIsotopicDistribution(targetMass,spectrum{1,jFile},"Window",6,"Tolerance",adjustedTolerance);
                    end
                end
                FeatureStruct(iFeature).isotopePattern = isotopePattern;
                FeatureStruct(iFeature).chargeState = chargeState;
            end
        end

        function featureStruct = calculateFormulaFromMass(obj,featureStruct)
            % finds chemical formula from monoisotopic mass
            MINSCORE = 600;
            tolerance = obj.betweenFileMassTolerance;
            tolUnit = obj.betweenFileMassUnit;

            parfor iFeature = 1:numel(featureStruct)
                currentMass = featureStruct(iFeature).mass_corrected;
                % determine Cl,Br and S counts from isotope distribution
                elementHits = detectIsotopicElements(featureStruct(iFeature).isotopePattern(:,1), featureStruct(iFeature).isotopePattern(:,2), featureStruct(iFeature).chargeState);
                maxCounts = [ceil(currentMass/12),...
                    floor(currentMass/1),...
                    floor(currentMass/79),...
                    floor(currentMass/35),...
                    floor(currentMass/19),...
                    floor(currentMass/127),...
                    floor(currentMass/14),...
                    floor(currentMass/16),...
                    floor(currentMass/31),...
                    floor(currentMass/32)];

                maxCounts(3) = elementHits.Br.count;
                maxCounts(4) = elementHits.Cl.count;
                maxCounts(10) = elementHits.S.count;

                % calculate formula
                if strcmp(tolUnit,"ppm")
                    toleranceAdjusted = tolerance * currentMass *10^-6;
                else
                    toleranceAdjusted = tolerance;
                end

                decomposition = fastMassDecomposition(currentMass,toleranceAdjusted,maxCounts);

                %filter nonsensical decompositions
                [evaluation,decomposition] = evaluateDecompositions(decomposition, currentMass);
                if isempty(decomposition)
                    continue
                end
                % compare feature isotope distribution to formula
                % distribution
                scores = findBestMatchingFormulaByIsotopeDistribution(decomposition,featureStruct(iFeature).isotopePattern,toleranceAdjusted);
                evaluation(scores < MINSCORE,:) = [];
                scores(scores < MINSCORE,:) = [];
                if ~isempty(evaluation)
                    [~,order] = sort(scores,"descend");
                    evaluation = evaluation(order,:);
                    featureStruct(iFeature).formulaEvaluation = evaluation;
                    featureStruct(iFeature).formula = evaluation.Formula(1);
                end
            end
        end

        %% Data handling
        function obj = initializeTemporaryFile(obj)
            %build TempDataFile
            obj.TempDataFile = tempname +".mat";
            obj.TempDataFileObj = matfile(obj.TempDataFile,Writable=true);
            %predefine Variables in .mat file
            obj.TempDataFileObj.ROICells  = {[]};
            obj.TempDataFileObj.TimeCells  = {[]};
            obj.TempDataFileObj.ROIMat = [];
            obj.TempDataFileObj.ROIMatProcessBLK  = {[]};
            obj.TempDataFileObj.ROIMatSystemBLK  = [];
            obj.TempDataFileObj.ROImzVec = [];
            obj.TempDataFileObj.timeVec  = [];
            % load MS1 Data from RawFile

            % MSData = load(obj.RawDataFile,"MSData");
            % MSData = MSData.MSData;
            MSData = obj.RawDataFileObj;
            nFiles = numel(MSData);
            timeCells = cell(nFiles,1);
            spectraCells = timeCells;
            for iFile = 1:nFiles
                timeCells{iFile,1} = [MSData(iFile).spectraMS1.rt]';
                spectraCells{iFile,1} = {MSData(iFile).spectraMS1.centroidedScan}';
            end
            obj.TempDataFileObj.ROICells = spectraCells;
            obj.TempDataFileObj.TimeCells = timeCells;
        end

        function obj = setOptimizationOptions(obj,optimizeMode,bayesOptions) %%%% WIP %%%%

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

    %% static methods block
    methods (Static)

        function IntegrationStruct = calculatePeakEntropy(IntegrationStruct)
            % Calculates Peak entropy for all peaks
            parfor iFeature = 1:length(IntegrationStruct)
                %check for no peaks, then skip iteration
                if ~isempty(IntegrationStruct(iFeature).peakStartLocation)
                    derivativ = diff(IntegrationStruct(iFeature).XIC(:,2));
                    probability = zeros(size(IntegrationStruct(iFeature).peakLocation));
                    for iPeak = 1:numel(probability)
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

        function integrationStruct = integratePeaks(integrationStruct,times)
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
                    %gather peak area
                    areas(jPeak,1) = trapz(eic(integrationStruct(iFeature).peakStartLocation(jPeak,1):integrationStruct(iFeature).peakEndLocation(jPeak,1)));
                end
                integrationStruct(iFeature).peakArea = areas;
                integrationStruct(iFeature).peakRetentionTime = times(integrationStruct(iFeature).peakLocation);
            end
        end

        function FeatureStruct = trimExtractedIonChromatograms(FeatureStruct)
            %% Trims the stored XIC to the location of the corresponding peak
            %
            nFile = width(FeatureStruct(1).peakLocations);
            allXic = {FeatureStruct.XIC}';
            allBorders = {FeatureStruct.peakBorders}';
            %preallocation
            trimedXIC = cell(length(FeatureStruct),nFile);

            parfor iFeat = 1:length(FeatureStruct)
                currentXIC = allXic{iFeat};
                currentBorders = allBorders{iFeat};
                for jFile = 1:nFile
                    %check for empty peak
                    if ~isnan(currentBorders(1,jFile))
                        trimedXIC{iFeat,jFile} = full(currentXIC(currentBorders(1,jFile):currentBorders(2,jFile),:));
                    else
                        continue
                    end
                end
            end
            %store
            for iFeat = 1:length(FeatureStruct)
                FeatureStruct(iFeat).XIC = trimedXIC(iFeat,:);
            end
        end

        function FeatureStruct = confirmSameFeatureByMS2(FeatureStruct)
            %% calculates the composit score within each feature
            % features with scores < 650 are split into a new feature

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
                spectra = denoiseScans(spectra,"threshold",0.025);

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

                    %check number of Peaks to remove
                    if isempty(idtoKeep) %all, keep first entry remove the rest
                        idtoSplit(1) = [];
                    elseif isempty(idtoSplit) %none, because of overlap
                        continue
                    else %some
                        id = any(idtoSplit == idtoKeep,1);
                        idtoSplit(id) = [];
                    end
                    splitFeatures = splitFeature(FeatureStruct(iFeat),idtoSplit);
                    %append
                    FeatureStruct(iFeat) = splitFeatures(1);
                    FeatureStruct = [FeatureStruct;splitFeatures(2:end)];
                end
            end
        end

        function FeatureStruct = confirmSameFeatureByIsotopeDistribution(FeatureStruct)
            %% compares isotopic pattern of each feature entry
            % features with bad scores are split into different features

            MINSCORE = 600;
            for iFeature = 1:length(FeatureStruct)
                isotopePattern = FeatureStruct(iFeature).isotopePattern;
                if sum(~cellfun("isempty",isotopePattern)) <= 1 %only one file with peak or no pattern
                    continue
                end

                %build index to original file
                originalFileID = 1:numel(isotopePattern);
                originalFileID(cellfun(@isempty,isotopePattern)) = [];
                isotopePattern(cellfun(@isempty,isotopePattern)) = [];

                alingedIsotopePattern = alignSpectra(isotopePattern,"normal","high","true");
                compoundScores = scoresWithinSet(alingedIsotopePattern);
                % rebuild original file
                for file = 1:numel(originalFileID)
                    compoundScores(compoundScores==file) = originalFileID(file);
                end

                if all(compoundScores(:,1) >= MINSCORE)
                    continue
                else %split feature
                    idtoKeep = compoundScores(compoundScores(:,1) >= MINSCORE,2:3);
                    idtoKeep = unique(idtoKeep);
                    idtoSplit = compoundScores(compoundScores(:,1) < MINSCORE,2:3);
                    idtoSplit = unique(idtoSplit);
                    if ~isempty(idtoKeep) & ~isempty(idtoKeep)
                        id = any(idtoSplit == idtoKeep,1);
                        idtoSplit(id) = [];
                    end
                    %check number of Peaks to remove
                    if isempty(idtoSplit)       %none, because of overlap
                        continue
                    elseif isempty(idtoKeep)    %all, keep first entry remove the rest

                        idtoSplit(1) = [];
                    end

                    splitFeatures = splitFeature(FeatureStruct(iFeature),idtoSplit);
                    %append
                    FeatureStruct(iFeature) = splitFeatures(1);
                    FeatureStruct = [FeatureStruct;splitFeatures(2:end)];
                end
            end
        end

        function FeatureStruct = averageIsotopePattern(FeatureStruct)
            %%  builds the average Isotopic pattern of each feature by aligning all patterns from each file and averaging them
            parfor iFeature = 1:height(FeatureStruct)
                isotopePattern = FeatureStruct(iFeature).isotopePattern;
                chargeState = FeatureStruct(iFeature).chargeState;
                %remove empty cells
                chargeState(cellfun(@isempty,isotopePattern)) = [];
                isotopePattern(cellfun(@isempty,isotopePattern)) = [];
                if ~isempty(chargeState) %average pattern and
                    FeatureStruct(iFeature).isotopePattern = alignSpectra(isotopePattern,"average","high","true");
                    FeatureStruct(iFeature).chargeState = median(chargeState,"all"); % only pattern of the same charge state, confirmation in previous step (confirmSameFeatureByIsotopeDistribution)
                else
                    FeatureStruct(iFeature).isotopePattern = [0,0];
                    FeatureStruct(iFeature).chargeState = 0;
                end
            end
            % filter features without isotope pattern
            FeatureStruct([FeatureStruct.chargeState] == 0) = [];
        end

        function FeatureStruct = correctMassByChargeState(FeatureStruct)

            function outputStruct = finalizeEISpectra(inputStruct)
                %% merges all found EI fragment spectra (all files) into an average spectrum
                outputStruct = inputStruct;
                for n = 1:length(inputStruct)
                    outputStruct(n).spectrumMS2 = alignSpectra(inputStruct(n).spectrumMS2,"average","low","true");
                end
            end

            %% uses charge state information to adjust the mass of each feature
            parfor iFeature = 1:height(FeatureStruct)
                if FeatureStruct(iFeature).chargeState > 0
                    FeatureStruct(iFeature).mass_corrected = FeatureStruct(iFeature).mass_measured * FeatureStruct(iFeature).chargeState;
                end
            end
        end
    end
end