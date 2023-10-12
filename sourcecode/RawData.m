classdef RawData
    % Class for storing group settings and performing funktions from Raw
    % data until Feature data stage
    properties
        %% Processing Parameters
        GroupName (1,1) string
        FileNames  (:,1) string
        Files (:,1) string {mustBeFile}
        BlankFiles (:,1) string {mustBeFile}
        % Main Processing Options
        MSPolarity (1,1) string {mustBeMember(MSPolarity,["positive","negative"])} = "positive"
        BLKSubtraction (1,1) logical = false
        Smoothing (1,1) logical = false
        BaseCorr (1,1) logical = false
        IsotopeFilter (1,1) logical = false
        AdductFilter (1,1) logical = false
        ContaminantFilter (1,1) logical = false
        ISTDCorr (1,1) logical = false
        MSalign (1,1) logical = false
        Peakalign (1,1) logical = false
        ScalingCorr (1,1) logical = false
        RTTol (1,1) double {mustBeFinite} = 5
        mzTol (1,1) double {mustBeFinite} = 0.001
        mzTolUnit (1,1) string {mustBeMember(mzTolUnit,["Da","ppm"])} = "Da"
        % ROIparameter
        thresh (1,1) double {mustBeInteger,mustBePositive} = 5000
        mzerror (1,1) double {mustBePositive} = 0.01
        mzErrorUnit (1,1) string {mustBeMember(mzErrorUnit,["Da","ppm"])} = "Da"
        minroi (1,1) double {mustBeInteger,mustBePositive} = 20
        Start (1,1) double {mustBeFinite} = 0
        End (1,1) double {mustBeFinite} = 1
        % Baseline Correction Parameters
        WindowSize (1,1) double {mustBeInteger,mustBePositive} = 200
        StepSize (1,1) double {mustBeInteger,mustBePositive} = 200
        RegressionMethod (1,1) string {mustBeMember(RegressionMethod,["pchip","linear","spline"])} = "pchip"
        EstimationMethod (1,1) string {mustBeMember(EstimationMethod,["quantile","em"])} = "em"
        SmoothMethod (1,1) string {mustBeMember(SmoothMethod,["none","lowess","loess","rlowess","rloess"])} = "none"
        QuantilVal (1,1) double {mustBeInRange(QuantilVal,0,1)} = 0.1
        %Golay Parameters
        FrameSize (1,1) double {mustBeInteger,mustBePositive} = 20
        Degree (1,1) double {mustBeInteger,mustBePositive} = 2
        % Internal Standard Data
        numISTD (1,1) double {mustBeInteger,mustBePositive} = 1
        ISDat (:,3) double
        MassCal (1,1) logical = false
        ISApply (1,1) string {mustBeMember(ISApply,["S&B","SOnly"])} = "SOnly"
        ISOrder (1,1) string {mustBeMember(ISOrder,["BlankIS","ISBlank"])} = "ISBlank"
        %Adduct Parameters
        CosSim (1,1) double {mustBeInRange(CosSim,0,1)} = 0.85
        AddSelectedPos (30,1) logical = false        %Structure: 1:12 Single Charged, 13:18 Dimers, 19:26 DoubleCharged, 27:30 TripleCharged
        AddSelectedNeg (16,1) logical = false       %Structure: 1:10 SingleCharged, 11:14 Dimers, 15 DoubleCharged, 16 TripleCharged
        NeutralSelectedSmol (18,1) logical = false
        NeutralSelectedCon (17,1) logical = false
        % mzalign Parameters
        mzEstimMethod (1,1) string {mustBeMember(mzEstimMethod,["histogram","regression"])} = "regression"
        mzCorrectionMethod (1,1) string {mustBeMember(mzCorrectionMethod,["nearest-neighbor","shortest-path"])} = "nearest-neighbor"
        mzQuantil (1,1) double {mustBeInRange(mzQuantil,0,1)} = 0.95
        % Peak Align Parameters
        maxshiftneg (1,1) double {mustBeInteger,mustBeNegative} = -50
        maxshiftpos (1,1) double {mustBeInteger,mustBePositive} = 50
        PulseWidth (1,1) double {mustBeInteger,mustBePositive} = 10
        WindowSizeRatio (1,1) double {mustBePositive} = 2.5
        SearchSpace (1,1) string {mustBeMember(SearchSpace,["regular","latin"])} = "regular"
        Iterations (1,1) double {mustBeInteger,mustBePositive} = 5
        GridSteps (1,1) double {mustBeInteger,mustBePositive} = 20
        % Scaling
        GroupScale (1,1) double {mustBePositive} = 1
        SampScale (:,1) double {mustBePositive} = 1
        % Integration and Filter
        EntropyFilter (1,1) logical = false
        minOccurence (1,1) double {mustBeInRange(minOccurence,0,1)} = 0.5
        minSignalNoise (1,1) double {mustBePositive} = 3
        maxWidth (1,1) double {mustBePositive} = 45
        minWidth (1,1) double {mustBePositive} = 2
        EvaluationParameter (1,1) string {mustBeMember(EvaluationParameter,["Area","Height"])} = "Area"
        %% DataStorage
        %Raw Data
        PreviewTICs (:,1) cell
        PreviewBPCs (:,1) cell
        PeakDataMS1 (:,1) cell
        PeakDataMSn (:,1) cell
        TimeDataMS1 (:,1) cell
        TimeDataMSn (:,1) cell
        Precursor   (:,1) cell
        CollisionType   (:,1) cell
        CollisionEnergy   (:,1) cell
        % processing temporaries
        ROICells    (:,1) cell
        TimeCells   (:,1) cell
        nScans      (:,1) double {mustBeInteger,mustBePositive}
        nScansPadded(:,1) double {mustBeInteger,mustBePositive}
        %FileInfos
        DataInfo    (:,4) double
        ScanFrequency (1,1) double
        %ROI Data
        ROIMat      (:,:) double
        ROIMatBLK   (:,:) double
        ROImzVec    (1,:) double
        timeVec     (:,1) double
        %IS Data
        ISValue     (:,:) double
        ISRT        (:,:) double
        ISMass      (1,:) double
        ISdelta     (1,:) double
        ISRTRange   (:,1) cell
        ISMZRange   (:,:) double
        % Number of removed Features
        SNFiltered       (1,1) double
        EntropyFiltered  (1,1) double
        MaxWidthFiltered (1,1) double
        MinWidthFiltered (1,1) double
        IsotopeFiltered  (1,1) double
        AdductFiltered   (1,1) double
        OccurenceFiltered(1,1) double
        MedianEntropy    (1,1) double
        EntropyStorage   (:,:) double
        Signal2NoiseStorage        (:,:) double
    end

    methods
        function obj = RawData(GroupCounter)
            %Construct an instance of this class
            obj.GroupName = "Group " + GroupCounter;
        end
        function [obj,polarity] = DataCheck(obj)
            % Check Data, number of Scans, Start/End Times
            %Check minimum number of scans
            FileLoc=[obj.Files;obj.BlankFiles];
            %remove empty
            idx=cellfun(@isempty,FileLoc);
            FileLoc(idx)=[];
            info=zeros(length(FileLoc),4);
            RetentionTimes = cell(length(FileLoc),1);
            TIC = cell(length(FileLoc),1);
            BPC = cell(length(FileLoc),1);
            FileInfo =  struct('Polarity','N/A',...
                'NumberOfScansMS1',[],...
                'NumberOfScansMSn',[],...
                'StartTime',[],...
                'EndTime',[]);
            parfor n=1:length(FileLoc)
                %check filetype
                test=strsplit(FileLoc(n),'.');
                test=test(end);
                switch test
                    case "mzML"
                        [FileInfo(n),RetentionTimes{n},TIC{n},BPC{n}] = mzMLinfo(FileLoc{n});
                    case "mzXML"
                        [FileInfo(n),RetentionTimes{n},TIC{n},BPC{n}] = mzXMLinfo(FileLoc{n});
                end
            end
            obj.PreviewTICs = TIC;
            obj.PreviewBPCs = BPC;
            obj.TimeDataMS1 = RetentionTimes;
            info(:,1)=[FileInfo.NumberOfScansMS1]+[FileInfo.NumberOfScansMSn];
            info(:,2)=[FileInfo.StartTime];
            info(:,3)=[FileInfo.EndTime];
            info(:,4)=(info(:,3)-info(:,2))./info(:,1);
            obj.Start=round(min(info(:,2)),1);
            obj.End=round(max(info(:,3)),1);
            obj.DataInfo = info;
            %check ms polarity of files
            pol = [FileInfo.Polarity];
            if all(pol=="positive")
                polarity = "positive";
                obj.MSPolarity = "positive";
            elseif all(pol=="negative")
                polarity = "negative";
                obj.MSPolarity = "negative";
            else
                polarity = "varied";
            end
        end

        function obj = ReadData(obj,DataLoc,Level)
            nFiles = size(DataLoc,1);
            %preallocation
            Peaks=cell(nFiles,1);
            times=cell(nFiles,1);
            if Level > 1
                PrecursorMass=cell(nFiles,1);
                CollisionForce=cell(nFiles,1);
                FragMethod=cell(nFiles,1);
            end
            parfor n=1:nFiles
                %filetype check
                FileType=strsplit(DataLoc(n),'.');
                FileType=FileType(end);
                switch FileType
                    case "mzML"
                        [Peaks{n,1},times{n,1},PrecursorMass{n,1},CollisionForce{n,1},FragMethod{n,1}] = readmzML(DataLoc{n},MSLevel=Level);
                    case "mzXML"
                        [Peaks{n,1},times{n,1},PrecursorMass{n,1},CollisionForce{n,1},FragMethod{n,1}] = readmzXML(DataLoc{n},MSLevel=Level);
                end
            end
            % compress peaklist by removing masses with intensity < 100;
            [Peaks,times] = DataCleanUp(Peaks,times);
            % store data
            if Level == 1
                obj.TimeDataMS1 = times;
                obj.PeakDataMS1 = Peaks;
            else
                obj.PeakDataMSn = Peaks;
                obj.TimeDataMSn = times;
                obj.Precursor = PrecursorMass;
                obj.CollisionEnergy = CollisionForce;
                obj.CollisionType = FragMethod;
            end
        end

        %% Data Processing
        function [obj,Output]=BatchProcess(obj)
            nFiles = size(obj.Files,1);
            nBLK = size(obj.BlankFiles,1);
            nData = nFiles+nBLK;
            FileLocs=[obj.Files;obj.BlankFiles];
            %remove possible empty
            id=cellfun(@isempty,FileLocs);
            FileLocs(id)=[];
            %check if files already loaded then skip loading stage
            if isempty(obj.PeakDataMS1) == true || size([obj.Files;obj.BlankFiles],1) ~= size(obj.PeakDataMS1,1)
                obj = ReadData(obj,FileLocs,1);
            end
            %remove scans outside RT range
            obj = obj.CutScansToSize;
            obj.nScans = cellfun(@numel,obj.TimeCells);
            if obj.MSalign == true
                obj = obj.AlignScans;
            end
            % ROI Search
            obj = obj.AutoROI;

            % Average BLK
            if obj.BLKSubtraction == true && nBLK > 1
                obj = obj.AverageBLK(nBLK);
                nData = size(obj.ROICells,1); % update number of matrices
            end
            if obj.ContaminantFilter == true
                obj = obj.removeContaminants;
            end

            % Baseline Corrention
            if obj.BaseCorr == true
                obj = obj.CorrectBaseline;
            end
            % Smoothing
            if obj.Smoothing == true
                obj = obj.SmoothPeaks;
            end

            % Peak Align
            if obj.Peakalign == true && nData > 1
                obj = obj.AlignPeaks;
            end

            if obj.BLKSubtraction == true % Separate Blank data from Sample data
                obj.ROIMatBLK=sparse(obj.ROICells{end});
                obj.ROICells(end)=[];
                obj.TimeCells(end)=[];
            end

            % subtract blank before IS normalization
            if obj.BLKSubtraction == true && obj.ISOrder == "BlankIS"
                for id=1:size(obj.ROICells,1)
                    obj.ROICells{id}=obj.ROICells{id}-obj.ROIMatBLK;
                    % set possible negative values to 0
                    obj.ROICells{id} = max(obj.ROICells{id},0);
                end
            end
            % pad arrays with Maximum peak width*3 Scans to eliminate
            % integration interference between matrices
            obj = obj.FinalizeROI;
            
            clearvars -except obj
            %% Integration Stage
            % Find and Integrate IS separate
            if obj.ISTDCorr == true
                obj = obj.IntegrateIS;
                obj = obj.findISTimeRanges;
                obj = obj.ISNormalize;
            end
            % BLK Subtraction after IS Correction
            if obj.BLKSubtraction == true && obj.ISOrder == "ISBlank"
                MSroi = mat2cell(obj.ROIMat,obj.nScansPadded);
                for id=1:size(MSroi,1)
                    MSroi{id}=MSroi{id}-padarray(obj.ROIMatBLK,size(MSroi{id},1)-size(obj.ROIMatBLK,1),0,'post');
                end
                obj.ROIMat = vertcat(MSroi{:});
                obj.ROIMat = max(obj.ROIMat,0);
                id = all(obj.ROIMat < obj.thresh,1);
                obj.ROIMat(:,id) = [];
                obj.ROImzVec(:,id) = [];
            end
            % mass correction
            if obj.MassCal == true
                obj = obj.findMZRanges;
                obj = obj.ISMassCorrection;
            end

            % Integrate all Peaks
            IDX = true(1,size(obj.ROIMat,2));
            IntegrationData = obj.CWTIntegrate(IDX);
            %Calculate number of removed features
            obj.MinWidthFiltered = sum(vertcat(IntegrationData{5,:}),"all");
            obj.MaxWidthFiltered = sum(vertcat(IntegrationData{6,:}),"all");
            obj.SNFiltered = sum(vertcat(IntegrationData{7,:}),"all");
            IntegrationData(5:7,:) = [];
            % remove ROI masses with no found peaks
            empt=cellfun(@isempty,IntegrationData(4,:));
            IntegrationData(:,empt)=[];
            obj.ROImzVec(empt)=[];
            obj.ROIMat(:,empt)=[];
            % calculate median entropy
            mEntropy=vertcat(IntegrationData{4,:});
            mEntropy(:,2)=[];
            obj.MedianEntropy=median(mEntropy,'omitnan');
            % Apply Entropy filter
            if obj.EntropyFilter == true
                [IntegrationData,obj.EntropyFiltered,EmptyColumns] = obj.FilterbyEntropy(IntegrationData,obj.MedianEntropy);
                obj.ROIMat(:,EmptyColumns)=[];
                obj.ROImzVec(EmptyColumns)=[];
            end
            IntegrationData = obj.AssignRT2SampleFile(IntegrationData);

            % remove isotopes and adducts
            if obj.IsotopeFilter == true
                [IntegrationData,obj] = obj.FilterIsotopes(IntegrationData);
            end

            if obj.AdductFilter == true
                [IntegrationData,obj] = obj.FilterAdducts(IntegrationData);
            end

            % Build Storage Arrays and filter by number of occurences
            [obj,Output] = obj.BuildStorageArrays(IntegrationData);
            % transform mz values from M+H+ / M-H+ to M
            switch obj.MSPolarity
                case "positive"
                    Output.FeatIdentifiers(:,1) = round(Output.FeatIdentifiers(:,1) - 1.007825,5);
                    obj.ROImzVec = round(obj.ROImzVec - 1.007825,5);
                case "negative"
                    Output.FeatIdentifiers(:,1) = round(Output.FeatIdentifiers(:,1) + 1.007825,5);
                    obj.ROImzVec = round(obj.ROImzVec + 1.007825,5);
            end
            obj.EntropyStorage = Output.EntropyStorage;
            obj.Signal2NoiseStorage = Output.Signal2NoiseStorage;
            Output.FeatIdentifiers(:,2) = round(Output.FeatIdentifiers(:,2),1);
            Output = obj.GroupAndSampleScaling(Output);
            Output.DataSize = size(Output.IntensityStorage,1);
            Output.FoundInGroup = repmat(obj.GroupName,size(Output.FeatIdentifiers,1),1);
            Output.SampleNames = obj.FileNames;
        end

        %% helper functions
        function obj = AverageBLK(obj,numBLK)
            ROICell = obj.ROICells;
            TimeCell = obj.TimeCells;
            maxScan=max(obj.nScans);
            BLKFiles = vertcat(ROICell{end-numBLK+1:end});
            BLKFiles = reshape(BLKFiles,maxScan,size(BLKFiles,2),numBLK);
            BLKFiles = mean(BLKFiles,3);
            BLKTimes = horzcat(TimeCell{end-numBLK+1:end});
            %replace 0 with NaN then ignore NaN in median
            %calculation
            BLKTimes(BLKTimes==0)=NaN;
            BLKTimes = median(BLKTimes,2,"omitnan");
            BLKTimes(isnan(BLKTimes))=0;
            ROICell(end-numBLK+1:end)=[];
            ROICell{end+1}=BLKFiles;
            TimeCell(end-numBLK+1:end)=[];
            TimeCell{end+1}=BLKTimes;
            obj.ROICells = ROICell;
            obj.TimeCells = TimeCell;
        end

        function obj = removeContaminants(obj)
            % Remove Contaminant Masses load correct Contaminant Masslist
            switch obj.MSPolarity
                case "positive"
                    Contaminants = load("MassListData.mat","ContaminantsPos");
                    Contaminants = Contaminants.ContaminantsPos;
                case "negative"
                    Contaminants = load("MassListData.mat","ContaminantsNeg");
                    Contaminants = Contaminants.ContaminantsNeg;
            end
            %calculate possible Contaminants
            switch obj.mzTolUnit
                case "Da"
                    isContaminant=abs(Contaminants-obj.ROImzVec) <= obj.mzTol;
                case "ppm"
                    isContaminant=abs(Contaminants-obj.ROImzVec)./obj.ROImzVec*10^6 <= obj.mzTol;
            end
            isContaminant=any(isContaminant,1);
            % remove contaminant columns from ROi mz list and MSroi
            % matrices
            obj.ROImzVec(isContaminant)=[];
            for id = 1:size(obj.ROICells,1)
                obj.ROICells{id}(:,isContaminant)=[];
            end
        end

        function obj = IntegrateIS(obj)
            % identify IS Vectors
            obj.ISMass = obj.ISDat(:,1)';
            [~,ISid] = min(abs(obj.ROImzVec-obj.ISMass'),[],2);
            %store difference and check tolerance
            obj.ISdelta = obj.ROImzVec(ISid)-obj.ISMass;
            switch obj.mzTolUnit
                case "Da"
                    dif = abs(obj.ISdelta) >= obj.mzTol;
                case "ppm"
                    dif = abs(obj.ISdelta)/obj.ROImzVec(ISid)*10^6 >= obj.mzTol;
            end
            %remove IS outside tolerance and throw warning
            if all(dif)
                % if no IS mass found, throw warning and exit
                header="Skipping ISTD Normalization";
                message = ["Reason: no Internal standard mass was found in this group","Check the specified m/z or increase the mass tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
                return
            elseif any(dif)
                % if some IS mass is not found, throw warning and continue
                header="Skipping ISTD No. ";
                for i=1:size(obj.ISMass,2)
                    if dif(i)==0
                        header = header + i + " ";
                    end
                end
                message = ["Reason: Internal standard mass was not found in this group","Check the specified m/z or increase the mass tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
            end
            ISid(dif)=[];
            obj.ISdelta(dif)=[];
            obj.ISMass(dif)=[];
            % extract relevant columns and perform Peak Picking and
            % Integration
            ISIntegrationData = obj.CWTIntegrate(ISid);
            % assign Peaks to Sample
            ISIntegrationData = obj.AssignRT2SampleFile(ISIntegrationData);
            % filter multiple Peaks for one sample
            for id=1:size(ISid)
                if numel(ISIntegrationData{3,id}(:,3)) ~= numel(unique(ISIntegrationData{3,id}(:,3)))
                    [~,idUniques,idAll] = unique(ISIntegrationData{3,id}(:,3));
                    for n=1:numel(idUniques)
                        testCase = idUniques(n) == idAll;
                        if sum(testCase)>1 % extract Peakdata of multiples, extract highest one
                            DataTemp = cellfun(@(x) x(testCase,:),ISIntegrationData(1:4,id),'UniformOutput',false);
                            ISIntegrationData(1:4,id) = cellfun(@(x) x(~testCase,:),ISIntegrationData(1:4,id),'UniformOutput',false);
                            [~,maxPeakid] = max(DataTemp{1,id});
                            ISIntegrationData(1:4,id) = cellfun(@(x,y) vertcat(y,x(maxPeakid,:)),DataTemp,ISIntegrationData(1:4,id),'UniformOutput',false);
                            [~,order] = sort(ISIntegrationData{3,id}(:,3),"ascend");
                            ISIntegrationData(1:4,id) = cellfun(@(x) x(order,:),ISIntegrationData(1:4,id),'UniformOutput',false);
                        end
                    end
                end
            end
            % check if RT Range is Correct and Remove Peaks outside range
            for id=1:size(ISid)
                RT = ISIntegrationData{3,id}(:,2);
                idx = ismembertol(RT,obj.ISDat(id,2),obj.ISDat(id,3),"DataScale",1);
                ISIntegrationData(1:4,id) = cellfun(@(x) x(idx,:), ISIntegrationData(1:4,id),'UniformOutput',false);
                ISSampleID = vertcat(ISIntegrationData{3,id}(:,3));
                RT(~idx) = 0;
                ISSampleID(~idx) = 0;
                % store RT data in object
                obj.ISRT(ISSampleID,id) = RT;
                switch obj.EvaluationParameter
                    case "Area"
                        ISVal = vertcat(ISIntegrationData{2,id}(:,1));
                        ISVal(~idx) = 0;
                        obj.ISValue(ISSampleID,id)=ISVal;
                    case "Height"
                        ISVal = vertcat(ISIntegrationData{1,id});
                        ISVal(~idx) = 0;
                        obj.ISValue(ISSampleID,id)=ISVal;
                end
            end
            % remove IS when not found in all samples and throw warning
            Check = any(obj.ISValue == 0,1);
            if any(Check) == true
                obj.ISRT(:,Check)=[];
                obj.ISValue(:,Check)=[];
                obj.ISdelta(:,Check)=[];
                obj.ISMass(:,Check)=[];
            end
            if all(Check)
                % if no IS mass found, throw warning and exit
                header="Skipping ISTD Normalization";
                message = ["Reason: no peaks found within the time tolerance","Check the specified retention time or increase the time tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
                return
            elseif any(Check)
                % if some IS mass is not found, throw warning and continue
                header="Skipping ISTD No. ";
                for i=1:size(Check,2)
                    if Check(i) == 1
                        header = header + i + ", ";
                    end
                end
                message = ["Reason: For one or more samples no suitable peak within the time tolerance was found","Check the specified retention time or increase the time tolerance"];
                fig = uifigure;
                uialert(fig,message,header,'Icon','warning');
            end
        end
        function IntResults = CWTIntegrate(obj,Index)
            Mat = obj.ROIMat(:,Index);
            maxSN = obj.minSignalNoise;
            % prepare wavelet filterbank
            MinPWDataPoints=floor(obj.minWidth/obj.ScanFrequency);
            MaxPWDataPoints=ceil(obj.maxWidth/obj.ScanFrequency);
            times = obj.timeVec;
            ScanFreq = obj.ScanFrequency;
            minSec = obj.minWidth;
            maxSec = obj.maxWidth;
            % calculate EIC derivatives and store as sparse
            smoothed = smoothdata(Mat,"gaussian","omitnan","SmoothingFactor",0.1);
            Noise = std(Mat-smoothed);
            Diff2 = zeros(length(times),size(Mat,2));
            Diff2(1:end-2,:) = diff(smoothed,2);
            numEIC = size(Mat,2);
            IntResults=cell(8,numEIC); %preallocate output
            for id=1:numEIC
                peaks = AutoCWT(Diff2(:,id),smoothed(:,id),ScanFreq,minSec,maxSec);
                % Correct Peak Borders
                peaks = CWTBorderCorrection(peaks,Mat(:,id),smoothed(:,id));
                [peaks,tempStorage] = FilterPeaks(peaks,MinPWDataPoints,MaxPWDataPoints,maxSN,Noise(:,id),Mat(:,id));
                %% integrate and store results
                tempStorage=FinalizeIntegrationOutput(peaks,tempStorage,Mat(:,id),times);
                for n=1:8
                    IntResults{n,id}=tempStorage{n,1};
                end
            end
        end

        function IntegrationData = AssignRT2SampleFile(obj,IntegrationData)
            test=cumsum(obj.nScansPadded)';
            nFiles = length(obj.nScansPadded);
            parfor n=1:size(IntegrationData,2)
                val = IntegrationData{3,n}(:,1);
                val = val < test;
                val = sum(val,2)-1;
                IntegrationData{3,n}(:,3) = abs(val-nFiles);
            end
        end

        function obj = ISNormalize(obj)
            %Gather relevant matrices
            if obj.ISApply == "S&B"
                Data = mat2cell(obj.ROIMat,obj.nScansPadded);
                Data{end+1} = obj.ROIMatBLK;
            else
                Data = mat2cell(obj.ROIMat,obj.nScansPadded);
            end

            % Retention time dependent IS normalization
            IntValues = obj.ISValue;
            RTRange = obj.ISRTRange;
            parfor id = 1:size(Data,1)
                SampleROI = Data{id};
                for nIS = 1:size(IntValues,2)
                    ISVal = IntValues(id,nIS);
                    ISRange = RTRange{id}(:,nIS);
                    SampleROI(ISRange(1):ISRange(2),:) = SampleROI(ISRange(1):ISRange(2),:)./ISVal;
                end
                Data{id} = SampleROI;
            end
            %store corrected Matrices back into object
            if obj.ISApply == "S&B"
                obj.ROIMatBLK = Data{end};
                Data(end) = [];
            end
            obj.ROIMat = vertcat(Data{:});
        end

        function obj = findISTimeRanges(obj)
            ISData = obj.ISRT;
            %split ROI matrix and timeVec into samples
            times = obj.timeVec;
            times=mat2cell(times,obj.nScansPadded,1);
            nScan=numel(times{1});
            for s=1:size(ISData,1) % sample loop
                ISrt=ISData(s,:);
                %RTval
                [ISrt,idx] = sort(ISrt,'ascend'); % sort rt
                TimeRange = arrayfun(@(i) mean(ISrt(i:i+1)),1:1:length(ISrt)-1)';
                TimeRange = round(TimeRange);
                rtborderEnds=[TimeRange;nScan];
                rtborderStarts=[1;TimeRange+1];
                TimeRange=[rtborderStarts rtborderEnds]';
                %resort Ranges
                unsorted = 1:length(ISrt);
                newIndRT(idx) = unsorted;
                obj.ISRTRange{s}=TimeRange(:,newIndRT);
            end
        end

        function obj = findMZRanges(obj)
            ISmz = obj.ISMass;
            mzVec = obj.ROImzVec;
            nMZ=numel(mzVec);
            %mzvalues
            [ISmz,idx] = sort(ISmz,'ascend'); % sort mz column
            MZRange = arrayfun(@(i) mean(ISmz(i:i+1)),1:1:length(ISmz)-1)'; %mean between two elements
            MZRange=MZRange-mzVec; %find indices of MZ values with minimum distance
            [~,MZRange]=min(abs(MZRange),[],2);
            mzborderEnds=[MZRange;nMZ]; % build ranges as 2 column matrix
            mzborderStarts=[1;MZRange+1];
            MZRange=[mzborderStarts mzborderEnds]';
            %resort Ranges
            unsorted = 1:length(ISmz);
            newIndMZ(idx) = unsorted;
            obj.ISMZRange=MZRange(:,newIndMZ);
        end

        function obj = ISMassCorrection(obj)
            CorrectionVector = zeros(size(obj.ROImzVec));
            deltas = obj.ISdelta;
            Ranges = obj.ISMZRange;
            for nIS = 1:numel(deltas)
                CorrectionVector(Ranges(1,nIS):Ranges(2,nIS)) = deltas(nIS);
            end
            obj.ROImzVec = obj.ROImzVec-CorrectionVector;
        end

        function [valuesFiltered,obj] = FilterAdducts(obj,IntegrationResults)
            %% AdductFilterAlgo Filters Adduct Peaks from Internal AriumMS integration results
            %   Calculates possible non Adduct (Base) m/z for each
            %   extracted m/z, then finds matching masses in original list.
            %   Peak lists of Adduct and Base m/z are then compared to have
            %   peaks with the same RT (within 2 sec). Matching peaks
            %   shapes are compared using Cosine Similarity (>=0.85
            %   default)
            switch obj.MSPolarity
                case "positive"
                    Rules = load("MassListData.mat","AddPosRules","NLossRules");
                    Rules = [Rules.AddPosRules(obj.AddSelectedPos);Rules.NLossRules([obj.NeutralSelectedSmol;obj.NeutralSelectedCon])];
                case "negative"
                    Rules = load("MassListData.mat","AddNegRules","NLossRules");
                    Rules = [Rules.AddNegRules(obj.AddSelectedNeg);Rules.NLossRules([obj.NeutralSelectedSmol;obj.NeutralSelectedCon])];
            end
            %Preparation
            minCosSim = obj.CosSim;
            mzVec = obj.ROImzVec;
            EICMat = obj.ROIMat;
            RTs=cellfun(@(X) X(:,2),IntegrationResults(3,:),'UniformOutput',false); %extract Retentiontimes
            ranges=cellfun(@(X) X(:,2:3),IntegrationResults(2,:),'UniformOutput',false); %extract Peak ranges

            %build possible Adduct mass list
            PossibleBaseMZ=zeros(length(Rules),length(mzVec));
            parfor n=1:length(Rules)
                PossibleBaseMZ(n,:)=Rules{n}(mzVec);
            end

            %compare case masses to initial mz list
            switch obj.mzTolUnit
                case "ppm"
                    PossibleBaseMZindex = cell(size(PossibleBaseMZ));
                    for n=1:size(mzVec,2)
                        [~,PossibleBaseMZindex(:,n)]=ismembertol(PossibleBaseMZ(:,n),mzVec(n),obj.mzTol,'ByRows',1,'DataScale',mzVec(n)*10^-6,'OutputAllIndices',true);
                    end
                case "Da"
                    [~,PossibleBaseMZindex]=ismembertol(PossibleBaseMZ,mzVec,obj.mzTol,'DataScale',1,'OutputAllIndices',true);
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
            obj.ROIMat(:,empt) = [];
            obj.ROImzVec(empt) = [];
        end

        function [valuesFiltered,obj] = FilterIsotopes(obj,IntegrationResults)
            %% IsotopeFilterAlgo Filters Isotope Peaks from Internal AriumMS integration results
            %
            %   Calculates possible non Isotope (Base) m/z forch each
            %   extracted m/z, then finds matching masses in original list.
            %   Peak lists of Isotope and Base m/z are then compared to
            %   have peaks with the same RT (within 2 sec). Matching peaks
            %   shapes are compared using Cosine Similarity (>=0.85), and
            %   Base m/z intensity adjusted by the relative Isotope
            %   occurence must be within 10% of the Isotope intensity

            %Preparation
            minCosSim = obj.CosSim;
            IsotopeRules = load("MassListData.mat","IsotopeRules");
            IsotopeRules = IsotopeRules.IsotopeRules;
            mz = obj.ROImzVec;
            EICMat = obj.ROIMat;
            RTs=cellfun(@(X) X(:,2),IntegrationResults(3,:),'UniformOutput',false); %extract Retentiontimes
            ranges=cellfun(@(X) X(:,2:3),IntegrationResults(2,:),'UniformOutput',false); %extract Peak ranges

            PossibleBaseMZ=mz-IsotopeRules(:,1); %build possible Isotope mass list
            %find possile BaseMZ values and gather Column index
            switch obj.mzTolUnit
                case "ppm"
                    PossibleBaseMZindex = cell(size(IsotopeRules,1),size(mz,2));
                    for n=1:size(mz,2)
                        [~,PossibleBaseMZindex(:,n)]=ismembertol(PossibleBaseMZ(:,n),mz(n),obj.mzTol,'ByRows',1,'DataScale',mz(n)*10^-6,'OutputAllIndices',true);
                    end
                case "Da"
                    [~,PossibleBaseMZindex]=ismembertol(PossibleBaseMZ,mz,obj.mzTol,'DataScale',1,'OutputAllIndices',true);
            end

            %build natural isotope occurence factor list
            IsotopeOccurenceFactors=cell(size(PossibleBaseMZindex));
            parfor n=1:size(PossibleBaseMZindex,1)
                IsotopeOccurenceFactors(n,:)=cellfun(@(X) ones(numel(X),1)*IsotopeRules(n,2), PossibleBaseMZindex(n,:), 'UniformOutput', false);
            end
            %%Check For matching Retention time
            %transform cell array to array
            rows=max(sum(cellfun(@numel,PossibleBaseMZindex)));
            idxMat=zeros(rows,length(mz),2);
            for k=1:length(mz)
                val=vertcat(PossibleBaseMZindex{:,k});
                idxMat(1:length(val),k,1)=val;
                val=vertcat(IsotopeOccurenceFactors{:,k});
                idxMat(1:length(val),k,2)=val;
            end
            %find matching RT indices
            matchRT=cell(size(idxMat(:,:,1)));
            for n=1:size(idxMat,2)
                for k=1:size(idxMat,1)
                    if idxMat(k,n)~=0
                        %check RT
                        [~,matchRT{k,n}]=ismembertol(RTs{n},RTs{idxMat(k,n,1)},1,'DataScale',1);
                    end
                end
            end
            %remove lists with only zeros
            check=cellfun(@(X) all(X(:)==0), matchRT);
            matchRT(check)={[]};
            check=repmat(check,1,1,2);
            idxMat(check)=0;
            %Preallocate storage variable of identified Isotopes variables
            IsotopeIndexCell = cellfun(@(X) false(size(X(:))), matchRT, 'UniformOutput', false);
            ColumnIndices = 1:1:size(matchRT,2);
            notemptyColumns = any(~cellfun(@isempty, matchRT),1);
            for n=ColumnIndices(notemptyColumns) %Isotope mz loop
                IsotopeEIC = EICMat(:,n);
                IsotopeRanges=ranges{1,n};
                RowIndices=1:1:size(matchRT,1);
                notEmptyRows = any(~cellfun(@isempty, matchRT(:,n)),2);
                for k=RowIndices(notEmptyRows) % possible Base mz loop
                    BaseEIC=EICMat(:,idxMat(k,n,1));    %base m/z EIC
                    BaseRanges=ranges{1,idxMat(k,n,1)}; %Base m/z peak ranges
                    CheckList=matchRT{k,n};             %List of Peaks that match Isotope rt
                    factor=idxMat(k,n,2);               %relative Intensity of Isotopes
                    %preallocate variables
                    MatchIndices=1:1:size(CheckList,1);
                    notZeroMatch = CheckList ~=0;
                    for p=MatchIndices(notZeroMatch)   % loop over matching Peaks
                        RangeA=IsotopeRanges(p,:);
                        RangeB=BaseRanges(CheckList(p,1),:);
                        EICA=IsotopeEIC;
                        EICB=BaseEIC;
                        %cut vectors to size
                        minFullRange = min([RangeA;RangeB],[],"all");
                        maxFullRange = max([RangeA;RangeB],[],"all");
                        EICA = EICA(minFullRange:maxFullRange);
                        EICB = EICB(minFullRange:maxFullRange);
                        tolA= max(EICA)+max(EICA)*0.10;
                        tolB= max(EICA)-max(EICA)*0.10;
                        %desision Cosine Similarity and mainPeak Intensity
                        %adjusted with relative Isotope Occurence is inside
                        %+-10% of isotopeintensity
                        IsotopeIndexCell{k,n}(p)=sum(EICB.*EICA)/(sqrt(sum(EICB.^2))*sqrt(sum(EICA.^2)))>=minCosSim & max(EICB)*factor<=tolA & max(EICB)*factor>=tolB;
                    end
                end
            end
            %remove identified Isotopes
            [valuesFiltered,numRemoved]=obj.Removify(IntegrationResults,IsotopeIndexCell);
            obj.IsotopeFiltered = numRemoved;
            %remove Empty columns
            empt=cellfun(@isempty,valuesFiltered(2,:));
            valuesFiltered(:,empt)=[];
            obj.ROIMat(:,empt) = [];
            obj.ROImzVec(empt) = [];
        end

        function [obj,Output] = BuildStorageArrays(obj,IntegrationResults)
            % Gather Data
            nFiles = size(obj.Files,1);
            minDataPoints = ceil(nFiles*obj.minOccurence);
            RTAssign = cellfun(@(x) x(:,2:3),IntegrationResults(3,:),'UniformOutput',false);
            TimeTolerance = obj.RTTol;
            switch obj.EvaluationParameter
                case "Height"
                    HeightOrArea = 1;
                case "Area"
                    HeightOrArea = 2;
            end
            nPeaks = cellfun(@(x) size(x,1),RTAssign);
            % remove cells with less peaks than required minimum
            idx = nPeaks<minDataPoints;

            %sum number of removed peaks peaks
            Removed = sum(nPeaks(idx),"all");
            % remove found cells
            RTAssign(idx)=[];
            IntegrationResults(:,idx)=[];
            obj.ROImzVec(idx)=[];
            obj.ROIMat(:,idx)=[];

            mzVector = obj.ROImzVec;

            % preallocate Storage CellArrays
            IntStorage = cell(size(mzVector));
            FeatId = cell(size(mzVector));
            XIC =  cell(size(mzVector));
            RTStorage = cell(size(mzVector));
            BorderStorage = cell(size(mzVector));
            Entropy_Storage = cell(size(mzVector));
            SNStorage = cell(size(mzVector));

            for n=1:size(RTAssign,2)
                Intensities = IntegrationResults{HeightOrArea,n}(:,1);
                Times = RTAssign{n}(:,1);
                LowerBorders = IntegrationResults{2,n}(:,2);
                UpperBorders = IntegrationResults{2,n}(:,3);
                Entropy = IntegrationResults{4,n}(:,1);
                SN = IntegrationResults{4,n}(:,2);
                SampleIndex = RTAssign{n}(:,2);
                XICvec = {[IntegrationResults{5,n},obj.timeVec]};
                %find unique Retention Times
                [UniqueTimes,IndexToUnique]=uniquetol(Times,TimeTolerance,'DataScale',1,'OutputAllIndices',true);

                % preallocate storage Matrices
                AvgTimeVec = zeros(size(UniqueTimes));
                IntMat = zeros(length(UniqueTimes),nFiles);
                TimesMat = zeros(length(UniqueTimes),nFiles);
                LowerBordersMat = zeros(length(UniqueTimes),nFiles);
                UpperBordersMat = zeros(length(UniqueTimes),nFiles);
                EntropyMat = zeros(length(UniqueTimes),nFiles);
                SNMat = zeros(length(UniqueTimes),nFiles);
                mzValue = repmat(mzVector(n),length(UniqueTimes),1);
                XICrep = repmat(XICvec,length(UniqueTimes),1);
                % uniqueRT loop
                for numRTs = 1:length(UniqueTimes)
                    AvgTimeVec(numRTs) = mean(Times(IndexToUnique{numRTs}));
                    IntMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = Intensities(IndexToUnique{numRTs});
                    TimesMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = Times(IndexToUnique{numRTs});
                    LowerBordersMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = LowerBorders(IndexToUnique{numRTs});
                    UpperBordersMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = UpperBorders(IndexToUnique{numRTs});
                    EntropyMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = Entropy(IndexToUnique{numRTs});
                    SNMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = SN(IndexToUnique{numRTs});
                end
                IntStorage{n} = IntMat;
                FeatId{n} = [mzValue,AvgTimeVec];
                XIC{n} = XICrep;
                RTStorage{n} = TimesMat;
                BordersMat = cat(3,LowerBordersMat,UpperBordersMat);
                BordersMat=mat2cell(BordersMat,ones(1,numel(UniqueTimes)),ones(1,nFiles),2);
                BordersMat = cellfun(@(x) squeeze(x), BordersMat, 'UniformOutput', false);
                BorderStorage{n} = BordersMat;
                Entropy_Storage{n} = EntropyMat;
                SNStorage{n} = SNMat;
            end
            Output.GroupName = obj.GroupName;
            Output.FeatIdentifiers = vertcat(FeatId{:});
            Output.XIC = vertcat(XIC{:});
            Output.IntensityStorage = vertcat(IntStorage{:});
            Output.RetentionTimeStorage = vertcat(RTStorage{:});
            Output.PeakBorderStorage = vertcat(BorderStorage{:});
            Output.EntropyStorage = vertcat(Entropy_Storage{:});
            Output.Signal2NoiseStorage = vertcat(SNStorage{:});

            %occurenceFilter
            idx = sum(Output.IntensityStorage ~= 0,2)<minDataPoints;
            Output.IntensityStorage(idx,:) = [];
            Output.FeatIdentifiers(idx,:) = [];
            Output.XIC(idx,:) = [];
            Output.RetentionTimeStorage(idx,:) = [];
            Output.PeakBorderStorage(idx,:) = [];
            Output.EntropyStorage(idx,:) = [];
            Output.Signal2NoiseStorage(idx,:) = [];
            obj.OccurenceFiltered = Removed + sum(idx);
        end
        function Output = GroupAndSampleScaling(obj,Output)
            %GroupScale
            Output.IntensityStorage = Output.IntensityStorage/obj.GroupScale;
            %SampleScale
            Output.IntensityStorage = Output.IntensityStorage./obj.SampScale';
        end
        function obj = CutScansToSize(obj)
            StartTime = obj.Start;
            EndTime = obj.End;
            tempPeakData = obj.PeakDataMS1;
            tempTimeData = obj.TimeDataMS1;
            parfor n = size(tempTimeData,1)
                idx = tempTimeData{n,1} < StartTime | tempTimeData{n,1} > EndTime
                tempPeakData{n,1}(idx)=[];
                tempTimeData{n,1}(idx)=[];
            end
            obj.ROICells = tempPeakData;
            obj.TimeCells = tempTimeData;
        end
        function obj=AutoROI(obj)
            %%AutoROI Performs fully automated ROI search and augmentation.
            Peaklist = obj.ROICells;
            Timelist = obj.TimeCells;
            Intthresh = obj.thresh;
            minroiSize = obj.minroi;
            ErrorUnit = obj.mzErrorUnit;
            Masserror = obj.mzerror;

            %preallocate cell arrays
            mzlist = cell(length(Peaklist),1);
            MSroilist = cell(length(Peaklist),1);
            %ROI search for every Sample
            parfor d = 1 : length(Peaklist)
                P= Peaklist{d,1};
                T= Timelist{d,1};
                nrows=length(P);
                [mzlist{d,1},MSroilist{d,1},~]=ROIpeaks2(P,Intthresh,Masserror,ErrorUnit,minroiSize,nrows,T);
            end
            if length(mzlist) == 1 %Skip Augmentation if only one Sample
                MSroi_end=MSroilist{1,1};
                mzroi_end=mzlist{1,1};
                time_end=Timelist{1,1};
            else
                for i = 2:size(Peaklist,1)
                    [MSroilist{1,1},mzlist{1,1},Timelist{1,1}] = MSroiaug2(MSroilist{1,1},MSroilist{i,1},mzlist{1,1},mzlist{i,1},Masserror,ErrorUnit,Intthresh,Timelist{1,1},Timelist{i,1});
                end
                MSroi_end=MSroilist{1,1};
                mzroi_end=mzlist{1,1};
                time_end=Timelist{1,1};
            end
            MSroi_end=MSroi_end-obj.thresh; %subtract intensity threshold
            MSroi_end=max(MSroi_end,0); %set every negative intensity to 0

            %split and pad matrices
            MSroi_end = mat2cell(MSroi_end,obj.nScans);
            time_end = mat2cell(time_end,obj.nScans);
            

            maxScan=max(obj.nScans);
            ScanNumbers = obj.nScans;
            parfor id = 1:size(MSroi_end,1)
                MSroi_end{id}= padarray(MSroi_end{id},maxScan-ScanNumbers(id),0,'post');
                time_end{id}= padarray(time_end{id},maxScan-ScanNumbers(id),0,'post');
            end
            obj.ROICells = MSroi_end;
            obj.TimeCells = time_end;
            obj.ROImzVec = mzroi_end;
        end
        function obj = AlignScans(obj)
            mzQuan = obj.mzQuantil;
            mzEstim = obj.mzEstimMethod;
            mzCorr = obj.mzCorrectionMethod;
            PeakCells = obj.ROICells;
            parfor id = 1:size(PeakCells,1)
                % perform Spectral Alignment
                [~, PeakCells{id}]= mspalign(PeakCells{id},'Quantile',mzQuan,'EstimationMethod',mzEstim,'CorrectionMethod',mzCorr,'ShowEstimation',false);
            end
            obj.ROICells = PeakCells;
        end
        function obj = CorrectBaseline(obj)
            WSize = obj.WindowSize;
            SSize = obj.StepSize;
            RegMethod = obj.RegressionMethod;
            EstMethod =obj.EstimationMethod;
            SmooMethod = obj.SmoothMethod;
            Quan = obj.QuantilVal;
            MSroi = obj.ROICells;
            time=obj.TimeCells;
            parfor id = 1:size(MSroi,1)
                oldSize=size(MSroi{id});
                %depad Array
                MSroiTemp = MSroi{id};
                [MSroiTemp,timeTemp] = depadArrays(MSroiTemp,time{id});
                MSroiTemp = msbackadj(timeTemp,MSroiTemp,'WindowSize',WSize,'StepSize',SSize,'RegressionMethod',RegMethod,'EstimationMethod',EstMethod,'SmoothMethod',SmooMethod,'QuantileValue',Quan,'PreserveHeights',true);
                %remove negativ, NaN and repad Array
                MSroiTemp=max(MSroiTemp,0);
                MSroiTemp(isnan(MSroiTemp))=0;
                [MSroi{id},time{id}] = repadArrays(MSroiTemp,timeTemp,oldSize);
                % set possible negative values to 0
                MSroi{id} = max(MSroi{id},0);
            end
            obj.ROICells = MSroi;
        end
        function obj = SmoothPeaks(obj)
            Frame = obj.FrameSize;
            Deg = obj.Degree;
            MSroi = obj.ROICells;
            time = obj.TimeCells;
            parfor id = 1:size(MSroi,1)
                %depad Array
                oldSize=size(MSroi{id});
                MSroiTemp = MSroi{id};
                [MSroiTemp,timeTemp] = depadArrays(MSroiTemp,time{id});
                MSroiTemp = mssgolay(timeTemp,MSroiTemp,'Span',Frame,'Degree',Deg);
                %remove negativ, NaN and repad Array
                MSroiTemp=max(MSroiTemp,0);
                MSroiTemp(isnan(MSroiTemp))=0;
                [MSroi{id},time{id}] = repadArrays(MSroiTemp,timeTemp,oldSize);
                % set possible negative values to 0
                MSroi{id} = max(MSroi{id},0);
            end
            obj.ROICells = MSroi;
        end
        function obj = AlignPeaks(obj)
            MSroi = obj.ROICells;
            time = obj.TimeCells;
            maxScan = max(obj.nScans);
            % rearrange matrices
            [splitVar,~] = cellfun(@size,time);
            test= cellfun(@(x) sum(x~=0),time);
            [~,test] = max(test);
            TimeVec = time{test};
            MSroi=vertcat(MSroi{:});
            time=vertcat(time{:});
            [~,id]=max(MSroi);
            refTime=time(id);
            ShiftVal = [obj.maxshiftneg,obj.maxshiftpos];
            PW = obj.PulseWidth;
            WSR = obj.WindowSizeRatio;
            I = obj.Iterations;
            GS = obj.GridSteps;
            SS = obj.SearchSpace;
            parfor n=1:size(MSroi,2)
                ROI=reshape(MSroi(:,n),maxScan,[]);
                ROI=msalign(TimeVec,ROI,refTime(n),'MaxShift',ShiftVal,...
                    'WidthOfPulses',PW,'WindowSizeRatio',WSR,'Iterations',...
                    I,'GridSteps',GS,'SearchSpace',SS);
                ROI(isnan(ROI))=0; %remove possible NaN
                MStemp{1,n} = reshape(ROI,[],1);
            end
            MSroi = cell2mat(MStemp);
            obj.ROICells = mat2cell(MSroi,splitVar);
            obj.TimeCells = mat2cell(time,splitVar);
        end
        function obj = FinalizeROI(obj)
            time = obj.TimeCells;
            MSroi = obj.ROICells;
            PaddedSize = zeros(size(MSroi));
            timeTemp=horzcat(time{:});
            timeTemp(any(timeTemp==0,2),:)=[];
            obj.ScanFrequency=mean(diff(timeTemp),'all');
            MaxPW=round(obj.maxWidth*1.5/obj.ScanFrequency);
            parfor i = 1:size(MSroi,1)
                MStemp=MSroi{i};
                MStemp=max(MStemp,0);
                MStemp(isnan(MStemp))=0;
                MStemp=padarray(MStemp,MaxPW,0,'post');
                MSroi{i}=MStemp;
                time{i}=padarray(time{i},MaxPW,0,'post');
                PaddedSize(i) = length(time{i});
            end
            obj.nScansPadded = PaddedSize;
            obj.nScans(length(PaddedSize)+1:end) = [];
            MStemp = vertcat(MSroi{:});
            %remove empty columns
            id = all(MStemp < obj.thresh,1);
            MStemp(:,id) = [];
            MStemp=max(MStemp,0);
            if obj.BLKSubtraction == true
                obj.ROIMatBLK(:,id) = [];
            end
            obj.ROImzVec(id) = [];
            obj.ROIMat = sparse(MStemp);
            obj.timeVec = round(vertcat(time{:}),1);
            %remove temporaries
            obj.TimeCells = [];
            obj.ROICells = [];
        end
    end
    %%
    methods(Static)

        function [IntegrationData,SumRemoved,empt] = FilterbyEntropy(IntegrationData,MedianEntropy)
            %preallocate number of removed Features
            SumRemoved = zeros(1,size(IntegrationData,2));
            for n=1:size(IntegrationData,2)
                idx=IntegrationData{4,n}(:,1)>MedianEntropy;
                IntegrationData{1,n}(idx,:)=[];
                IntegrationData{2,n}(idx,:)=[];
                IntegrationData{3,n}(idx,:)=[];
                IntegrationData{4,n}(idx,:)=[];
                SumRemoved(n) = sum(idx);
            end
            % remove empty features and store
            empt=cellfun(@isempty,IntegrationData(2,:));
            IntegrationData(:,empt)=[];
            SumRemoved = sum(SumRemoved,"all");
        end

        function [values,sumRemoved]=Removify(values,found)
            %removes confirmed Adduct from Input list and stores number of
            %removed peaks
            sumRemoved = zeros(1,size(found,2));
            for i=1:size(found,2)
                idx = any(horzcat(found{:,i}),2);
                if isempty(idx)==false
                    values{1,i}(idx,:)=[];
                    values{2,i}(idx,:)=[];
                    values{3,i}(idx,:)=[];
                    values{4,i}(idx,:)=[];
                    sumRemoved(i)=sum(idx);
                else
                    sumRemoved(i) = 0;
                end
            end
            sumRemoved = sum(sumRemoved,"all");
        end

    end
end