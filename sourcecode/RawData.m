classdef RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        %% Processing Parameters
        FileNames           (:,1) string
        Files               (:,1) string
        BlankFiles          (:,1) string
        % Main Processing Options
        MSFileType          (1,1) string {mustBeMember(MSFileType,["profile","centroid"])} = "profile"
        BLKSubtraction      (1,1) logical = false
        Smoothing           (1,1) logical = false
        BaseCorr            (1,1) logical = false
        IsotopeFilter       (1,1) logical = false
        AdductFilter        (1,1) logical = false
        ContaminantFilter   (1,1) logical = false
        ISTDCorr            (1,1) logical = false
        MSalign             (1,1) logical = false
        Peakalign           (1,1) logical = false
        ScalingCorr         (1,1) logical = false

        % ROI parameter
        thresh              (1,1) double {mustBeInteger,mustBePositive} = 5000
        mzerror             (1,1) double {mustBePositive} = 0.01
        mzErrorUnit         (1,1) string {mustBeMember(mzErrorUnit,["Da","ppm"])} = "Da"
        minroi              (1,1) double {mustBeInteger,mustBePositive} = 20
        Start               (1,1) double {mustBeFinite} = 0
        End                 (1,1) double {mustBeFinite} = 1
        % Baseline Correction Parameters
        WindowSize          (1,1) double {mustBeFinite,mustBePositive} = 20
        StepSize            (1,1) double {mustBeFinite,mustBePositive} = 20
        RegressionMethod    (1,1) string {mustBeMember(RegressionMethod,["pchip","linear","spline"])} = "pchip"
        EstimationMethod    (1,1) string {mustBeMember(EstimationMethod,["quantile","em"])} = "em"
        SmoothMethod        (1,1) string {mustBeMember(SmoothMethod,["none","lowess","loess"])} = "none"
        QuantilVal          (1,1) double {mustBeInRange(QuantilVal,0,1)} = 0.1
        %Golay Parameters
        FrameSize           (1,1) double {mustBeInteger,mustBePositive} = 20
        Degree              (1,1) double {mustBeInteger,mustBePositive} = 2
        % Internal Standard Data
        numISTD             (1,1) double {mustBeInteger,mustBePositive} = 1
        ISDat               (:,3) double
        MassCal             (1,1) logical = false
        ISApply             (1,1) string {mustBeMember(ISApply,["S&B","SOnly"])} = "SOnly"
        ISOrder             (1,1) string {mustBeMember(ISOrder,["BlankIS","ISBlank"])} = "ISBlank"
        %Adduct Parameters
        CosSim              (1,1) double {mustBeInRange(CosSim,0,1)} = 0.85
        AddSelectedPos      (30,1) logical = false        %Structure: 1:12 Single Charged, 13:18 Dimers, 19:26 DoubleCharged, 27:30 TripleCharged
        AddSelectedNeg      (16,1) logical = false       %Structure: 1:10 SingleCharged, 11:14 Dimers, 15 DoubleCharged, 16 TripleCharged
        NeutralSelectedSmol (18,1) logical = false
        NeutralSelectedCon  (17,1) logical = false
        % mzalign Parameters
        mzEstimMethod       (1,1) string {mustBeMember(mzEstimMethod,["histogram","regression"])} = "regression"
        mzCorrectionMethod  (1,1) string {mustBeMember(mzCorrectionMethod,["nearest-neighbor","shortest-path"])} = "nearest-neighbor"
        mzQuantil           (1,1) double {mustBeInRange(mzQuantil,0,1)} = 0.99
        % Peak Align Parameters
        maxshiftneg         (1,1) double {mustBeFinite,mustBePositive} = 20
        maxshiftpos         (1,1) double {mustBeFinite,mustBePositive} = 20
        PulseWidth          (1,1) double {mustBeFinite,mustBePositive} = 2
        WindowSizeRatio     (1,1) double {mustBePositive} = 2.5
        SearchSpace         (1,1) string {mustBeMember(SearchSpace,["regular","latin"])} = "regular"
        Iterations          (1,1) double {mustBeInteger,mustBePositive} = 5
        GridSteps           (1,1) double {mustBeInteger,mustBePositive} = 20
        % Scaling
        GroupScale          (1,1) double {mustBePositive} = 1
        SampScale           (:,1) double {mustBePositive} = 1
        % Integration and Filter
        EvaluationParameter (1,1) string {mustBeMember(EvaluationParameter,["Area","Height"])} = "Area"
        minWidth            (1,1) double {mustBePositive} = 2
        maxWidth            (1,1) double {mustBePositive} = 45
        minSignalNoise      (1,1) double {mustBePositive} = 3
        minOccurence        (1,1) double {mustBeInRange(minOccurence,0,1)} = 0.5
        mzTol               (1,1) double {mustBeFinite} = 0.001
        mzTolUnit           (1,1) string {mustBeMember(mzTolUnit,["Da","ppm"])} = "Da"
        RTTol               (1,1) double {mustBeFinite} = 5
        entropyFilter       (1,1) logical = false
        entropyStrength     (1,1) string {mustBeMember(entropyStrength,["strict","medium","lax"])} = "medium"
        %% DataStorage
        RawDataFile         string
        RawDataFileObj      (1,1)
        TempDataFile        string
        TempDataFileObj     (1,1)
        ROIDataFile         string
        ROIDataFileObj      (1,1)
        % processing variables
        nScans              (:,1) double {mustBeInteger,mustBePositive}
        nScansPadded        (:,1) double {mustBeInteger,mustBePositive}
        %FileInfos
        DataInfo            (1,:) struct
        ScanFrequency       (1,1) double
        %IS Data
        ISValue             (:,:) double
        ISRT                (:,:) double
        ISMass              (1,:) double
        ISMassFound         (1,:) double
        ISdelta             (:,:) double
        mzCorrectionFcn     (1,1)
        % Number of removed Features
        SNFiltered          (1,1) double
        EntropyFiltered     (1,1) double
        MaxWidthFiltered    (1,1) double
        MinWidthFiltered    (1,1) double
        IsotopeFiltered     (1,1) double
        AdductFiltered      (1,1) double
        OccurenceFiltered   (1,1) double
        MedianEntropy       (1,1) double

        %testing variables
        Output
    end

    methods
        function obj = RawData
            %Construct an instance of this class
            obj.RawDataFile = tempname +".mat";
            obj.RawDataFileObj = matfile(obj.RawDataFile,Writable=true);

            %predefine Variables in .mat file
            obj.RawDataFileObj.PreviewTICs = {[]};
            obj.RawDataFileObj.PreviewBPCs = {[]};
            obj.RawDataFileObj.PreviewTimes = {[]};

            obj.RawDataFileObj.PeakDataMS1 = {[]};
            obj.RawDataFileObj.TimeDataMS1 = {[]};

            obj.RawDataFileObj.PeakDataMSn = {[]};
            obj.RawDataFileObj.TimeDataMSn = {[]};
            obj.RawDataFileObj.Precursor = {[]};
            obj.RawDataFileObj.CollisionEnergy = {[]};
            obj.RawDataFileObj.CollisionType = {[]};

            obj.ROIDataFile = tempname +".mat";
        end

        function obj = DataCheck(obj)
            % Check Data, number of Scans, Start/End Times
            %Check minimum number of scans
            FileLoc=[obj.Files;obj.BlankFiles];
            %remove empty
            idx=cellfun(@isempty,FileLoc);
            FileLoc(idx)=[];
            RetentionTimes = cell(length(FileLoc),1);
            TIC = cell(length(FileLoc),1);
            BPC = cell(length(FileLoc),1);
            polarityCells = cell(length(FileLoc),1);
            FileInfo =  struct('NumberOfScansMS1',[],...
                'NumberOfScansMSn',[],...
                'StartTime',[],...
                'EndTime',[]);
            parfor n=1:length(FileLoc)
                %check filetype
                test=strsplit(FileLoc(n),'.');
                test=test(end);
                switch test
                    case "mzML"
                        [FileInfo(n),RetentionTimes{n},TIC{n},BPC{n},polarityCells{n}] = mzMLinfo(FileLoc{n});
                    case "mzXML"
                        [FileInfo(n),RetentionTimes{n},TIC{n},BPC{n},polarityCells{n}] = mzXMLinfo(FileLoc{n});
                end
            end
            obj.RawDataFileObj.PreviewTICs = TIC;
            obj.RawDataFileObj.PreviewBPCs = BPC;
            obj.RawDataFileObj.PreviewTimes = RetentionTimes;
            obj.RawDataFileObj.polarity = polarityCells;
            obj.RawDataFileObj.PeakDataMS1 ={[]};

            % calculate Scan Frequency [Hz]
            scanFrq = [FileInfo.NumberOfScansMS1]./([FileInfo.EndTime]-[FileInfo.StartTime]);
            scanFrq = num2cell(scanFrq);
            [FileInfo.ScanFrequenceMS1] = scanFrq{:};
            scanFrq = [FileInfo.NumberOfScansMSn]./([FileInfo.EndTime]-[FileInfo.StartTime]);
            scanFrq = num2cell(scanFrq);
            [FileInfo.ScanFrequenceMSn] = scanFrq{:};
            %store data
            obj.DataInfo = FileInfo;
            obj.Start=round(min([FileInfo.StartTime]),1);
            obj.End=round(max([FileInfo.EndTime]),1);
        end

        function obj = ReadData(obj,DataLoc,Level)
            nFiles = size(DataLoc,1);
            %preallocation
            Peaks=cell(nFiles,1);
            times=cell(nFiles,1);
            PrecursorMass=cell(nFiles,1);
            CollisionForce=cell(nFiles,1);
            FragMethod=cell(nFiles,1);
            fileType = obj.MSFileType;
            %check if DataCheck was performed
            if ~isfield(obj.RawDataFileObj,"polarity")
                obj = obj.DataCheck;
            end
            
            polarities = obj.RawDataFileObj.polarity;
            parfor n=1:nFiles
                peakTemp = [];
                timeTemp = [];
                %filetype check
                FileType=strsplit(DataLoc(n),'.');
                FileType=FileType(end);
                switch FileType
                    case "mzML"
                        [peakTemp,timeTemp,PrecursorMass{n,1},CollisionForce{n,1},FragMethod{n,1}] = readmzML(DataLoc{n},MSLevel=Level);
                    case "mzXML"
                        [peakTemp,timeTemp,PrecursorMass{n,1},CollisionForce{n,1},FragMethod{n,1}] = readmzXML(DataLoc{n},MSLevel=Level);
                end
                
                % when profile data then centroid scans
                if fileType == "profile"
                    peakTemp = CentroidScans(peakTemp);
                else
                    [peakTemp,timeTemp] = DataCleanUp(peakTemp,timeTemp);
                end
                %convert from pseudo molecular mass to molecular mass
                peakTemp = ConvertScans2MolecularMass(peakTemp,polarities{n});
                Peaks{n,1} = peakTemp;
                times{n,1} = timeTemp;
            end
            if Level == 1
                obj.RawDataFileObj.TimeDataMS1 = times;
                obj.RawDataFileObj.PeakDataMS1 = Peaks;
            else
                %remove cells with no MSn data
                idx = cellfun(@isempty,Peaks);
                Peaks(idx,:) = [];
                times(idx,:) = [];
                PrecursorMass(idx,:) = [];
                CollisionForce(idx,:) = [];
                FragMethod(idx,:) = [];
                [obj.RawDataFileObj.PeakDataMSn,obj.RawDataFileObj.TimeDataMSn,obj.RawDataFileObj.Precursor,obj.RawDataFileObj.CollisionEnergy,obj.RawDataFileObj.CollisionType] = obj.MS2CleanUp(Peaks,times,PrecursorMass,CollisionForce,FragMethod);
            end
        end

        %% Data Processing
        function [Output,obj]=BatchProcess(obj,varargin)
            %check if old results exist and delete them
            if isfile(obj.ROIDataFile)
                delete(obj.ROIDataFile)
            end

            if numel(varargin) == 2
                mode = varargin{1};
                bayesOptions = varargin{2};
                obj = obj.SetOptimizationOptions(mode,bayesOptions);
            elseif isscalar(varargin)
                error("Wrong number of inputs")
            end

            nFiles = size(obj.Files,1);
            FileLocs = obj.Files;
            nData = nFiles;
            if obj.BLKSubtraction == true
                nBLK = size(obj.BlankFiles,1);
                nData = nFiles+nBLK;
                FileLocs=[FileLocs;obj.BlankFiles];
            end
            %remove possible empty
            id=cellfun(@isempty,FileLocs);
            FileLocs(id)=[];
            nData = nData-sum(id);

            %check if files already loaded then skip loading stage
            test = obj.RawDataFileObj.PeakDataMS1(1,1);
            if isempty(test{1,1}) || size([obj.Files;obj.BlankFiles],1) ~= size(obj.RawDataFileObj.PeakDataMS1,1)
                obj = obj.ReadData(FileLocs,1);
            end
            clearvars test FileLocs id

            %build TempDataFile
            obj.TempDataFile = tempname +".mat";
            obj.TempDataFileObj = matfile(obj.TempDataFile,Writable=true);
            %predefine Variables in .mat file
            obj.TempDataFileObj.ROICells  = {[]};
            obj.TempDataFileObj.TimeCells  = {[]};
            obj.TempDataFileObj.ROIMat = [];
            obj.TempDataFileObj.ROIMatBLK  = [];
            obj.TempDataFileObj.ROImzVec = [];
            obj.TempDataFileObj.timeVec  = [];

            %remove scans outside RT range
            obj = obj.CutScansToSize;
            
             % remove isotopes and adducts
            if obj.IsotopeFilter == true
                obj = obj.FilterIsotopesScanStage;
            end

            obj.nScans = cellfun(@numel,obj.TempDataFileObj.TimeCells);
            if obj.MSalign == true
                obj = obj.AlignScans("batch");
            end

            % ROI Search
            obj = obj.AutoROI("batch");

            % Average BLK
            if obj.BLKSubtraction == true && nBLK > 1
                obj = obj.AverageBLK(nBLK);
                nData = size(obj.TempDataFileObj.ROICells,1); % update number of matrices
            end
            if obj.ContaminantFilter == true
                obj = obj.removeContaminants;
            end

            % Baseline Correction
            if obj.BaseCorr == true
                obj = obj.CorrectBaseline("batch");
            end
            % Smoothing
            if obj.Smoothing == true
                obj = obj.SmoothPeaks("batch");
            end

            % Peak Align
            if obj.Peakalign == true && nData > 1
                obj = obj.AlignPeaks("batch");
            end

            if obj.BLKSubtraction == true % Separate Blank data from Sample data
                tempBLK = obj.TempDataFileObj.ROICells(end,1);
                obj.TempDataFileObj.ROIMatBLK=sparse(tempBLK{:});
                obj.TempDataFileObj.ROICells(end)=[];
                obj.TempDataFileObj.TimeCells(end)=[];
            end

            % subtract blank before IS normalization
            if obj.BLKSubtraction == true && obj.ISOrder == "BlankIS"
                peakCells = obj.TempDataFileObj.ROICells;
                BLKMat = obj.TempDataFileObj.ROIMatBLK;
                parfor id=1:size(peakCells,1)
                    peakCells{id,1}=peakCells{id,1}-BLKMat;
                    % set possible negative values to 0
                    peakCells{id,1} = max(peakCells{id,1},0);
                end
                obj.TempDataFileObj.ROICells = peakCells;
            end
            % pad arrays with Maximum peak width*3 Scans to eliminate
            % integration interference between matrices
            obj = obj.FinalizeROI;

            clearvars -except obj
            %% Integration Stage
            % Find and Integrate IS separate
            if obj.ISTDCorr == true
                obj = obj.IntegrateIS;
                if ~isempty(obj.ISValue)
                    obj = obj.ISNormalize;
                end
            end
            % BLK Subtraction after IS Correction
            if obj.BLKSubtraction == true && obj.ISOrder == "ISBlank"
                MSroi = mat2cell(obj.TempDataFileObj.ROIMat,obj.nScansPadded);
                MatBLK = obj.TempDataFileObj.ROIMatBLK;
                parfor id=1:size(MSroi,1)
                    MSroi{id,1}=MSroi{id,1}-padarray(MatBLK,size(MSroi{id,1},1)-size(MatBLK,1),0,'post');
                end
                MSroi = vertcat(MSroi{:});
                MSroi = max(MSroi,0);
                id = all(MSroi >= obj.thresh,1);
                obj.TempDataFileObj.ROIMat = MSroi(:,id);
                obj.TempDataFileObj.ROImzVec(:,~id) = [];
            end
            % mass correction
            if obj.MassCal == true && ~isempty(obj.ISValue)
                obj = obj.ISMassCorrection;
            end

            % Integrate all Peaks
            IDX = true(1,size(obj.TempDataFileObj.ROIMat,2));
            IntegrationData = obj.CWTIntegrate(IDX);
            %Calculate number of removed features
            obj.MinWidthFiltered = sum(vertcat(IntegrationData{5,:}),"all");
            obj.MaxWidthFiltered = sum(vertcat(IntegrationData{6,:}),"all");
            obj.SNFiltered = sum(vertcat(IntegrationData{7,:}),"all");
            IntegrationData(5:7,:) = [];
            % remove ROI masses with no found peaks
            IDX=cellfun(@isempty,IntegrationData(4,:));
            IntegrationData(:,IDX)=[];
            obj.TempDataFileObj.ROImzVec(IDX)=[];
            obj.TempDataFileObj.ROIMat(:,IDX)=[];
            % calculate median entropy
            mEntropy=vertcat(IntegrationData{4,:});
            if ~isempty(mEntropy)
                mEntropy(:,2)=[];
            end
            obj.MedianEntropy=median(mEntropy,'omitnan');
            % Apply Entropy filter
            if obj.entropyFilter == true
                [IntegrationData,obj.EntropyFiltered,EmptyColumns] = obj.FilterbyEntropy(IntegrationData,obj.MedianEntropy);
                obj.TempDataFileObj.ROIMat(:,EmptyColumns)=[];
                obj.TempDataFileObj.ROImzVec(EmptyColumns)=[];
            end
            IntegrationData = obj.AssignRT2SampleFile(IntegrationData);

            % remove adducts
            if obj.AdductFilter == true
                [IntegrationData,obj] = obj.FilterAdducts(IntegrationData);
            end

            % Build Storage Arrays and filter by number of occurences
            [Output,obj] = obj.BuildStorageArrays(IntegrationData);

            %check for empty Output
            if isempty(Output.FeatIdentifiers)
                Output.FeatIdentifiers(1,1:2) = 0;
            end
            Output.FeatIdentifiers(:,2) = round(Output.FeatIdentifiers(:,2),1);
            Output = obj.GroupAndSampleScaling(Output);
            Output.DataSize = size(Output.IntensityStorage,1);
            Output.FoundInGroup = repmat(obj.GroupName,size(Output.FeatIdentifiers,1),1);
            Output.SampleNames = obj.FileNames;
            obj.Output = Output;
            if ~exist("mode","var") %save results if batch mode
                %build ResultDataFile
                obj.ROIDataFileObj = matfile(obj.ROIDataFile,Writable=true);
                %store Results
                obj.ROIDataFileObj.ROIMat = obj.TempDataFileObj.ROIMat;
                obj.ROIDataFileObj.ROIMatBLK  = obj.TempDataFileObj.ROIMatBLK;
                obj.ROIDataFileObj.ROImzVec = obj.TempDataFileObj.ROImzVec;
                obj.ROIDataFileObj.timeVec  = obj.TempDataFileObj.timeVec;
            end
            %delete Temprorary file
            delete(obj.TempDataFile)
            obj.TempDataFile = "";
        end

        function obj = SetOptimizationOptions(obj,OptimizeMode,bayesOptions)

            switch OptimizeMode

                case "MainOptions"
                    if ismember("intThresh",bayesOptions.Properties.VariableNames)
                        obj.thresh = bayesOptions.intThresh;
                    end
                    if ismember("mzerror",bayesOptions.Properties.VariableNames)
                        obj.mzerror = bayesOptions.mzerror;
                    end
                    if ismember("minRoi",bayesOptions.Properties.VariableNames)
                        obj.minroi = bayesOptions.minRoi;
                    end
                    if ismember("minPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.minPeakWidth)
                        obj.minWidth = bayesOptions.minPeakWidth;
                    end
                    if ismember("maxPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.maxPeakWidth)
                        obj.maxWidth = bayesOptions.maxPeakWidth;
                    end
                    if ismember("minSN",bayesOptions.Properties.VariableNames)
                        obj.minSignalNoise = bayesOptions.minSN;
                    end
                    if ismember("mzTol",bayesOptions.Properties.VariableNames)
                        obj.mzTol = bayesOptions.mzTol;
                    end
                    if ismember("timeTol",bayesOptions.Properties.VariableNames)
                        obj.RTTol = bayesOptions.timeTol;
                    end
                    if ismember("entropyFilter",bayesOptions.Properties.VariableNames)
                        obj.entropyFilter = bayesOptions.entropyFilter == "true";
                    end
                    if ismember("blankCorrection",bayesOptions.Properties.VariableNames)
                        obj.BLKSubtraction = bayesOptions.blankCorrection == "true";
                    end
                    if ismember("contaminantFilter",bayesOptions.Properties.VariableNames)
                        obj.ContaminantFilter = bayesOptions.contaminantFilter == "true";
                    end
                    if ismember("isotopeFilter",bayesOptions.Properties.VariableNames)
                        obj.IsotopeFilter = bayesOptions.isotopeFilter == "true";
                    end
                    if ismember("MSAlign",bayesOptions.Properties.VariableNames)
                        obj.MSalign = bayesOptions.MSAlign == "true";
                    end
                    if ismember("peakAlignment",bayesOptions.Properties.VariableNames)
                        obj.Peakalign = bayesOptions.peakAlignment == "true";
                    end
                    if ismember("baselineCorrection",bayesOptions.Properties.VariableNames)
                        obj.BaseCorr = bayesOptions.baselineCorrection == "true";
                    end
                    if ismember("smoothing",bayesOptions.Properties.VariableNames)
                        obj.Smoothing = bayesOptions.smoothing == "true";
                    end

                case "SubParameters"
                    if ismember("entropyStrength",bayesOptions.Properties.VariableNames)
                        obj.entropyStrength = bayesOptions.entropyStrength;
                    end
                    %MSAlign parameters
                    if ismember("mzEstimMethod",bayesOptions.Properties.VariableNames)
                        obj.mzEstimMethod = bayesOptions.mzEstimMethod;
                    end
                    if ismember("mzCorrectionMethod",bayesOptions.Properties.VariableNames)
                        obj.mzCorrectionMethod = bayesOptions.mzCorrectionMethod;
                    end
                    if ismember("mzQuantil",bayesOptions.Properties.VariableNames)
                        obj.mzQuantil = bayesOptions.mzQuantil;
                    end
                    %PeakAlign parameters
                    if ismember("maxShiftneg",bayesOptions.Properties.VariableNames)
                        obj.maxshiftneg = bayesOptions.maxShiftneg;
                    end
                    if ismember("maxShiftpos",bayesOptions.Properties.VariableNames)
                        obj.maxshiftpos = bayesOptions.maxShiftpos;
                    end
                    if ismember("pulseWidth",bayesOptions.Properties.VariableNames)
                        obj.PulseWidth = bayesOptions.pulseWidth;
                    end
                    if ismember("iterations",bayesOptions.Properties.VariableNames)
                        obj.Iterations = bayesOptions.iterations;
                    end
                    if ismember("searchSpace",bayesOptions.Properties.VariableNames)
                        obj.SearchSpace = bayesOptions.searchSpace;
                    end
                    if ismember("gridSteps",bayesOptions.Properties.VariableNames)
                        obj.GridSteps = bayesOptions.gridSteps;
                    end
                    %baseline parameters
                    %always set smoothing to none
                    obj.SmoothMethod = "none";
                    if ismember("windowSize",bayesOptions.Properties.VariableNames)
                        obj.WindowSize = bayesOptions.windowSize;
                    end
                    if ismember("stepSize",bayesOptions.Properties.VariableNames)
                        obj.StepSize = bayesOptions.stepSize;
                    end
                    if ismember("regressionMethod",bayesOptions.Properties.VariableNames)
                        obj.RegressionMethod = bayesOptions.regressionMethod;
                    end
                    if ismember("estimationMethod",bayesOptions.Properties.VariableNames)
                        obj.EstimationMethod = bayesOptions.estimationMethod;
                    end
                    if ismember("quantile",bayesOptions.Properties.VariableNames)
                        obj.QuantilVal = bayesOptions.quantile;
                    end
                    %Smoothing parameters
                    if ismember("frameSize",bayesOptions.Properties.VariableNames)
                        obj.FrameSize = bayesOptions.frameSize;
                    end
                    if ismember("polyDegree",bayesOptions.Properties.VariableNames)
                        obj.Degree = bayesOptions.polyDegree;
                    end


                case {"Full","Custom"}
                    if ismember("intThresh",bayesOptions.Properties.VariableNames)
                        obj.thresh = bayesOptions.intThresh;
                    end
                    if ismember("mzerror",bayesOptions.Properties.VariableNames)
                        obj.mzerror = bayesOptions.mzerror;
                    end
                    if ismember("minRoi",bayesOptions.Properties.VariableNames)
                        obj.minroi = bayesOptions.minRoi;
                    end
                    if ismember("minPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.minPeakWidth)
                        obj.minWidth = bayesOptions.minPeakWidth;
                    end
                    if ismember("maxPeakWidth",bayesOptions.Properties.VariableNames) && ~isnan(bayesOptions.maxPeakWidth)
                        obj.maxWidth = bayesOptions.maxPeakWidth;
                    end
                    if ismember("minSN",bayesOptions.Properties.VariableNames)
                        obj.minSignalNoise = bayesOptions.minSN;
                    end
                    if ismember("mzTol",bayesOptions.Properties.VariableNames)
                        obj.mzTol = bayesOptions.mzTol;
                    end
                    if ismember("timeTol",bayesOptions.Properties.VariableNames)
                        obj.RTTol = bayesOptions.timeTol;
                    end
                    if ismember("entropyFilter",bayesOptions.Properties.VariableNames)
                        obj.entropyFilter = bayesOptions.entropyFilter == "true";
                        if bayesOptions.entropyFilter == "true"
                            if ismember("entropyStrength",bayesOptions.Properties.VariableNames)
                                obj.entropyStrength = bayesOptions.entropyStrength;
                            end
                        end
                    end
                    if ismember("blankCorrection",bayesOptions.Properties.VariableNames)
                        obj.BLKSubtraction = bayesOptions.blankCorrection == "true";
                    end
                    if ismember("contaminantFilter",bayesOptions.Properties.VariableNames)
                        obj.ContaminantFilter = bayesOptions.contaminantFilter == "true";
                    end
                    if ismember("isotopeFilter",bayesOptions.Properties.VariableNames)
                        obj.IsotopeFilter = bayesOptions.isotopeFilter == "true";
                    end
                    if ismember("MSAlign",bayesOptions.Properties.VariableNames)
                        obj.MSalign = bayesOptions.MSAlign == "true";
                    end
                    if ismember("peakAlignment",bayesOptions.Properties.VariableNames)
                        obj.Peakalign = bayesOptions.peakAlignment == "true";
                    end
                    if ismember("baselineCorrection",bayesOptions.Properties.VariableNames)
                        obj.BaseCorr = bayesOptions.baselineCorrection == "true";
                    end
                    if ismember("smoothing",bayesOptions.Properties.VariableNames)
                        obj.Smoothing = bayesOptions.smoothing == "true";
                    end
                    %MS Alignment
                    if ismember("MSAlign",bayesOptions.Properties.VariableNames)
                        obj.MSalign = bayesOptions.MSAlign == "true";
                        if bayesOptions.MSAlign == "true"
                            if ismember("mzEstimMethod",bayesOptions.Properties.VariableNames)
                                obj.mzEstimMethod = bayesOptions.mzEstimMethod;
                            end
                            if ismember("mzCorrectionMethod",bayesOptions.Properties.VariableNames)
                                obj.mzCorrectionMethod = bayesOptions.mzCorrectionMethod;
                            end
                            if ismember("mzQuantil",bayesOptions.Properties.VariableNames)
                                obj.mzQuantil = bayesOptions.mzQuantil;
                            end
                        end
                    end
                    if ismember("peakAlignment",bayesOptions.Properties.VariableNames)
                        obj.Peakalign = bayesOptions.peakAlignment == "true";
                        if bayesOptions.peakAlignment == "true"
                            %Peak Alignment
                            if ismember("maxShiftneg",bayesOptions.Properties.VariableNames)
                                obj.maxshiftneg = bayesOptions.maxShiftneg;
                            end
                            if ismember("maxShiftpos",bayesOptions.Properties.VariableNames)
                                obj.maxshiftpos = bayesOptions.maxShiftpos;
                            end
                            if ismember("pulseWidth",bayesOptions.Properties.VariableNames)
                                obj.PulseWidth = bayesOptions.pulseWidth;
                            end
                            if ismember("iterations",bayesOptions.Properties.VariableNames)
                                obj.Iterations = bayesOptions.iterations;
                            end
                            if ismember("searchSpace",bayesOptions.Properties.VariableNames)
                                obj.SearchSpace = bayesOptions.searchSpace;
                            end
                            if ismember("gridSteps",bayesOptions.Properties.VariableNames)
                                obj.GridSteps = bayesOptions.gridSteps;
                            end
                        end
                    end
                    %Baseline Correction
                    if ismember("baselineCorrection",bayesOptions.Properties.VariableNames)
                        obj.BaseCorr = bayesOptions.baselineCorrection == "true";
                        if bayesOptions.baselineCorrection == "true"
                            %always set smoothing to none
                            obj.SmoothMethod = "none";
                            if ismember("windowSize",bayesOptions.Properties.VariableNames)
                                obj.WindowSize = bayesOptions.windowSize;
                            end
                            if ismember("stepSize",bayesOptions.Properties.VariableNames)
                                obj.StepSize = bayesOptions.stepSize;
                            end
                            if ismember("regressionMethod",bayesOptions.Properties.VariableNames)
                                obj.RegressionMethod = bayesOptions.regressionMethod;
                            end
                            if ismember("estimationMethod",bayesOptions.Properties.VariableNames)
                                obj.EstimationMethod = bayesOptions.estimationMethod;
                            end
                            if ismember("smoothingMethod",bayesOptions.Properties.VariableNames)
                                obj.SmoothMethod = bayesOptions.smoothingMethod;
                            end
                            if ismember("quantile",bayesOptions.Properties.VariableNames)
                                obj.QuantilVal = bayesOptions.quantile;
                            end
                        end
                    end
                    % Smoothing
                    if ismember("smoothing",bayesOptions.Properties.VariableNames)
                        obj.Smoothing = bayesOptions.smoothing == "true";
                        if bayesOptions.smoothing == "true"
                            if ismember("frameSize",bayesOptions.Properties.VariableNames)
                                obj.FrameSize = bayesOptions.frameSize;
                            end
                            if ismember("polyDegree",bayesOptions.Properties.VariableNames)
                                obj.Degree = bayesOptions.polyDegree;
                            end
                        end
                    end
            end
        end

        %% helper functions
        function obj = AverageBLK(obj,numBLK)
            ROICell = obj.TempDataFileObj.ROICells;
            TimeCell = obj.TempDataFileObj.TimeCells;
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
            obj.TempDataFileObj.ROICells = ROICell;
            obj.TempDataFileObj.TimeCells = TimeCell;
        end

        function obj = removeContaminants(obj)
            % Remove Contaminant Masses load correct Contaminant Masslist
            polarity = obj.RawDataFileObj.polarity;
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
                    Contaminants = load("MassListData.mat","ContaminantsPos");
                    Contaminants = Contaminants.ContaminantsPos;
                case "negative"
                    Contaminants = load("MassListData.mat","ContaminantsNeg");
                    Contaminants = Contaminants.ContaminantsNeg;
                case "both"
                    Contaminants = load("MassListData.mat","ContaminantsPos","ContaminantsNeg");
                    Contaminants = unique([Contaminants.ContaminantsPos;Contaminants.ContaminantsNeg]);
            end
            %calculate possible Contaminants
            switch obj.mzTolUnit
                case "Da"
                    isContaminant=abs(Contaminants-obj.TempDataFileObj.ROImzVec) <= obj.mzTol;
                case "ppm"
                    isContaminant=abs(Contaminants-obj.TempDataFileObj.ROImzVec)./obj.TempDataFileObj.ROImzVec*10^6 <= obj.mzTol;
            end
            isContaminant=any(isContaminant,1);
            % remove contaminant columns from ROi mz list and MSroi
            % matrices
            obj.TempDataFileObj.ROImzVec(isContaminant)=[];
            tempCell = obj.TempDataFileObj.ROICells;
            for id = 1:size(obj.TempDataFileObj.ROICells,1)
                tempCell{id,1}(:,isContaminant)=[];
            end
            obj.TempDataFileObj.ROICells = tempCell;
        end

        function obj = IntegrateIS(obj)
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
            ISIntegrationData = obj.CWTIntegrate(ISid);
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
            ISData = obj.BuildStorageArrays(ISIntegrationData,obj.ISMassFound);

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

        function IntResults = CWTIntegrate(obj,varargin)
            if isscalar(varargin)
                Index = varargin{1};
                Mat = obj.TempDataFileObj.ROIMat;
                Mat = Mat(:,Index);
            else
                Mat = sum(obj.TempDataFileObj.ROIMat,2);
            end

            minSN = obj.minSignalNoise;

            % prepare wavelet filter-bank
            MinPWDataPoints=floor(obj.minWidth/obj.ScanFrequency);
            MaxPWDataPoints=ceil(obj.maxWidth/obj.ScanFrequency);
            times = obj.TempDataFileObj.timeVec;

            % calculate EIC derivatives and store as sparse
            smoothed = smoothdata(Mat,"gaussian","omitnan","SmoothingFactor",0.1);
            Noise = std(Mat-smoothed);
            Diff2 = zeros(length(times),size(Mat,2));
            Diff2(1:end-2,:) = diff(smoothed,2);
            numEIC = size(Mat,2);
            IntResults=cell(8,numEIC); %preallocate output
            FilterBank = cwtfilterbank("SignalLength",size(Diff2,1),"WaveletParameters",[3 4],"VoicesPerOctave",8,"SamplingPeriod",seconds(obj.ScanFrequency),"PeriodLimits",[seconds(obj.minWidth) seconds(obj.maxWidth)]);% prepare wavelet filterbank
            parfor id=1:numEIC
                peaks = AutoCWT(Diff2(:,id),smoothed(:,id),FilterBank);
                % Correct Peak Borders
                peaks = CWTBorderCorrection(peaks,Mat(:,id),smoothed(:,id));
                [peaks,tempStorage] = FilterPeaks(peaks,MinPWDataPoints,MaxPWDataPoints,minSN,Noise(:,id),Mat(:,id));
                % integrate and store results
                IntResults(:,id) = FinalizeIntegrationOutput(peaks,tempStorage,Mat(:,id),times);
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
                Data = mat2cell(obj.TempDataFileObj.ROIMat,obj.nScansPadded);
                Data{end+1} = obj.TempDataFileObj.ROIMatBLK;
            else
                Data = mat2cell(obj.TempDataFileObj.ROIMat,obj.nScansPadded);
            end

            timeVectors = mat2cell(obj.TempDataFileObj.timeVec,obj.nScansPadded);
            
            % normalize Intensities
            if isscalar(obj.ISValue) %only one IS
                value = obj.ISValue;
                parfor id = 1:size(Data,1)
                    Data{id} = Data{id}./value;
                end

            else % multiple IS - Retention time dependent IS normalization
                retentionTimes = obj.ISRT;
                intensities = obj.ISValue;
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
            if obj.ISApply == "S&B"
                obj.TempDataFileObj.ROIMatBLK = Data{end};
                Data(end) = [];
            end
            obj.TempDataFileObj.ROIMat = vertcat(Data{:});
        end

        function obj = ISMassCorrection(obj)
            if isscalar(obj.ISdelta) % single IS constant correction
                obj.TempDataFileObj.ROImzVec = obj.TempDataFileObj.ROImzVec-obj.ISdelta;

            else  % multiple IS - m/z dependent correction
                %% fit correction function
                [xData, yData] = prepareCurveData(obj.ISMass, obj.ISdelta);
                % Set up fittype and options.
                ft = 'pchipinterp';
                opts = fitoptions( 'Method', 'PchipInterpolant' );
                opts.ExtrapolationMethod = 'nearest';
                % Fit model to data.
                obj.mzCorrectionFcn = fit(xData,yData,ft,opts);

                % build mass correction vector and subtract from ROI masses
                mzCorrectionVector = obj.mzCorrectionFcn(obj.TempDataFileObj.ROImzVec);
                obj.TempDataFileObj.ROImzVec = obj.TempDataFileObj.ROImzVec - mzCorrectionVector;
            end
        end

        function [valuesFiltered,obj] = FilterAdducts(obj,IntegrationResults)
            %% AdductFilterAlgo Filters Adduct Peaks from Internal AriumMS integration results
            %   Calculates possible non Adduct (Base) m/z for each
            %   extracted m/z, then finds matching masses in original list.
            %   Peak lists of Adduct and Base m/z are then compared to have
            %   peaks with the same RT (within 2 sec). Matching peaks
            %   shapes are compared using Cosine Similarity (>=0.85
            %   default)

            mzTolVal = obj.mzTol;

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
                    Rules = [Rules.AddPosRules(obj.AddSelectedPos);Rules.NLossRules([obj.NeutralSelectedSmol;obj.NeutralSelectedCon])];
                case "negative"
                    Rules = load("MassListData.mat","AddNegRules","NLossRules");
                    Rules = [Rules.AddNegRules(obj.AddSelectedNeg);Rules.NLossRules([obj.NeutralSelectedSmol;obj.NeutralSelectedCon])];
                case "both"
                    Rules = load("MassListData.mat","AddNegRules","AddNegRules","NLossRules");
                    Rules = [Rules.AddPosRules(obj.AddSelectedPos);Rules.AddNegRules(obj.AddSelectedNeg);Rules.NLossRules([obj.NeutralSelectedSmol;obj.NeutralSelectedCon])];
            end
            %Preparation
            minCosSim = obj.CosSim;
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
            switch obj.mzTolUnit
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
        function obj = FilterIsotopesScanStage(obj)
            tolerance = obj.mzerror;
            tolUnit = obj.mzErrorUnit;
            tempPeakData = obj.TempDataFileObj.ROICells;
            for n = 1:size(tempPeakData,1)
                tempPeakData{n,1} = InScanIsotopeFilter(tempPeakData{n,1},tolerance,tolUnit);
            end
            obj.TempDataFileObj.ROICells = tempPeakData;
        end
        
        function [Output,obj] = BuildStorageArrays(obj,IntegrationResults,varargin)
            TimeTolerance = obj.RTTol;
            nFiles = size(obj.nScans,1);
            if isscalar(varargin)
                mzVector = varargin{1};
                minDataPoints = nFiles;
                isISIntegration = true;
            else
                mzVector = obj.TempDataFileObj.ROImzVec;
                minOcc = obj.minOccurence;
                minDataPoints = ceil(nFiles*minOcc);
                isISIntegration = false;
            end

            % Gather Data

            RTAssign = cellfun(@(x) x(:,2:3),IntegrationResults(3,:),'UniformOutput',false);
            switch obj.EvaluationParameter
                case "Height"
                    HeightOrArea = 1;
                case "Area"
                    HeightOrArea = 2;
            end
            nPeaks = cellfun(@(x) size(x,1),RTAssign);
            % remove cells with less peaks than required minimum
            idx = nPeaks<minDataPoints;

            %sum number of removed peaks
            Removed = sum(nPeaks(idx),"all");
            % remove cells with fever then required peaks
            RTAssign(idx)=[];
            IntegrationResults(:,idx)=[];
            mzVector(:,idx)=[];
            if isISIntegration == false
                obj.TempDataFileObj.ROImzVec(:,idx) = [];
                obj.TempDataFileObj.ROIMat(:,idx) = [];
            end

            % preallocate Storage CellArrays
            IntStorage = cell(size(mzVector));
            FeatId = cell(size(mzVector));
            XIC =  cell(size(mzVector));
            RTStorage = cell(size(mzVector));
            Entropy_Storage = cell(size(mzVector));
            SNStorage = cell(size(mzVector));
            timeVec = obj.TempDataFileObj.timeVec;

            intensities = IntegrationResults(HeightOrArea,:);
            lowerBorders = IntegrationResults(2,:);
            upperBorders = IntegrationResults(2,:);
            entropyAndSN = IntegrationResults(4,:);
            XICvec = IntegrationResults(5,:);

            parfor n=1:size(RTAssign,2)
                localIntensity = intensities{n}(:,1)
                localLowerBorder = lowerBorders{n}(:,2);
                localUpperBorder = upperBorders{n}(:,3);
                localEntropy = entropyAndSN{n}(:,1);
                localSN = entropyAndSN{n}(:,2);
                localXIC = [XICvec{n},timeVec];

                Times = RTAssign{n}(:,1);
                SampleIndex = RTAssign{n}(:,2);
                
                %find unique Retention Times
                [UniqueTimes,IndexToUnique] = uniquetol(Times,TimeTolerance,'DataScale',1,'OutputAllIndices',true);

                % preallocate storage Matrices
                AvgTimeVec = zeros(size(UniqueTimes));
                IntMat = zeros(length(UniqueTimes),nFiles);
                TimesMat = zeros(length(UniqueTimes),nFiles);
                LowerBordersMat = zeros(length(UniqueTimes),nFiles);
                UpperBordersMat = zeros(length(UniqueTimes),nFiles);
                EntropyMat = zeros(length(UniqueTimes),nFiles);
                SNMat = zeros(length(UniqueTimes),nFiles);
                mzValue = repmat(mzVector(n),length(UniqueTimes),1);
                % uniqueRT loop

                for numRTs = 1:length(UniqueTimes)
                    AvgTimeVec(numRTs) = mean(Times(IndexToUnique{numRTs}));
                    IntMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = localIntensity(IndexToUnique{numRTs});
                    TimesMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = Times(IndexToUnique{numRTs});
                    LowerBordersMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = localLowerBorder(IndexToUnique{numRTs});
                    UpperBordersMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = localUpperBorder(IndexToUnique{numRTs});
                    EntropyMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = localEntropy(IndexToUnique{numRTs});
                    SNMat(numRTs,SampleIndex(IndexToUnique{numRTs})) = localSN(IndexToUnique{numRTs});
                end
                IntStorage{n} = IntMat;
                FeatId{n} = [mzValue,AvgTimeVec];
                RTStorage{n} = TimesMat;
                BordersMat = cat(3,LowerBordersMat,UpperBordersMat);
                BordersMat=mat2cell(BordersMat,ones(1,numel(UniqueTimes)),ones(1,nFiles),2);
                BordersMat = cellfun(@(x) squeeze(x), BordersMat, 'UniformOutput', false);
                XIC{n} = ExtractXIC(localXIC,BordersMat);
                Entropy_Storage{n} = EntropyMat;
                SNStorage{n} = SNMat;
            end
            Output.GroupName = obj.GroupName;
            Output.FeatIdentifiers = vertcat(FeatId{:});
            Output.XIC = vertcat(XIC{:});
            Output.IntensityStorage = vertcat(IntStorage{:});
            Output.RetentionTimeStorage = vertcat(RTStorage{:});
            Output.EntropyStorage = vertcat(Entropy_Storage{:});
            Output.Signal2NoiseStorage = vertcat(SNStorage{:});
            % duplicate row filter
            [Output.FeatIdentifiers,idx] = unique(Output.FeatIdentifiers,'rows','stable');
            Output.IntensityStorage = Output.IntensityStorage(idx,:);
            Output.XIC = Output.XIC(idx,:);
            Output.RetentionTimeStorage = Output.RetentionTimeStorage(idx,:);
            Output.EntropyStorage = Output.EntropyStorage(idx,:);
            Output.Signal2NoiseStorage = Output.Signal2NoiseStorage(idx,:);

            %occurenceFilter
            idx = sum(Output.IntensityStorage ~= 0,2)<minDataPoints;
            Output.IntensityStorage(idx,:) = [];
            Output.FeatIdentifiers(idx,:) = [];
            Output.XIC(idx,:) = [];
            Output.RetentionTimeStorage(idx,:) = [];
            Output.EntropyStorage(idx,:) = [];
            Output.Signal2NoiseStorage(idx,:) = [];
             if isISIntegration == false
                obj.OccurenceFiltered = Removed + sum(idx);
            end
            
        end

        function Output = GroupAndSampleScaling(obj,Output)
            if ~isempty(Output.IntensityStorage)
                %GroupScale
                Output.IntensityStorage = Output.IntensityStorage/obj.GroupScale;
                %SampleScale
                Output.IntensityStorage = Output.IntensityStorage./obj.SampScale';
            end
        end

        function obj = CutScansToSize(obj)
            StartTime = obj.Start;
            EndTime = obj.End;
            tempPeakData = obj.RawDataFileObj.PeakDataMS1;
            tempTimeData = obj.RawDataFileObj.TimeDataMS1;

            parfor n = 1:size(tempPeakData,1)
                idx = tempTimeData{n,1} < StartTime | tempTimeData{n,1} > EndTime;
                tempPeakData{n,1}(idx)=[];
                tempTimeData{n,1}(idx)=[];
            end
            obj.TempDataFileObj.ROICells = tempPeakData;
            obj.TempDataFileObj.TimeCells = tempTimeData;
        end

        function obj = AutoROI(obj,modeFlag)
            %%AutoROI Performs fully automated ROI search and augmentation.
            switch modeFlag
                case "preview"
                    peakList = obj.TempDataFileObj.ROICells(1,1);
                    timeList = obj.TempDataFileObj.TimeCells(1,1);
                case "batch"
                    peakList = obj.TempDataFileObj.ROICells;
                    timeList = obj.TempDataFileObj.TimeCells;
            end
            intThresh = obj.thresh;
            minroiSize = obj.minroi;
            errorUnit = obj.mzErrorUnit;
            massError = obj.mzerror;

            %preallocate cell arrays
            mzlist = cell(length(peakList),1);
            MSroilist = cell(length(peakList),1);
            %ROI search for every Sample
            parfor d = 1 : length(peakList)
                P= peakList{d,1};
                T= timeList{d,1};
                [mzlist{d,1},MSroilist{d,1},~]=ROIpeaks3(P,intThresh,massError,errorUnit,minroiSize,T);
            end
            if isscalar(mzlist) %Skip Augmentation if only one Sample
                MSroi_end=MSroilist{1,1};
                mzroi_end=mzlist{1,1};
                time_end=timeList{1,1};
            else
                for i = 2:size(peakList,1)
                    [MSroilist{1,1},mzlist{1,1},timeList{1,1}] = MSroiaug3(MSroilist{1,1},MSroilist{i,1},mzlist{1,1},mzlist{i,1},massError,errorUnit,intThresh,timeList{1,1},timeList{i,1});
                end
                MSroi_end=MSroilist{1,1};
                mzroi_end=mzlist{1,1};
                time_end=timeList{1,1};
            end
            MSroi_end=MSroi_end-obj.thresh; %subtract intensity threshold
            MSroi_end=max(MSroi_end,0); %set every negative intensity to 0

            if modeFlag ~= "preview"
                %split and pad matrices
                outROI = mat2cell(MSroi_end,obj.nScans);
                outTime = mat2cell(time_end,obj.nScans);
                maxScan=max(obj.nScans);
                ScanNumbers = obj.nScans;
                parfor id = 1:size(outROI,1)
                    outROI{id,1} = padarray(outROI{id,1},maxScan-ScanNumbers(id,1),0,'post');
                    outTime{id,1} = padarray(outTime{id,1},maxScan-ScanNumbers(id,1),0,'post');
                end
            else
                outROI{1,1} = MSroi_end;
                outTime{1,1} = time_end;
            end
            obj.TempDataFileObj.ROICells = outROI;
            obj.TempDataFileObj.TimeCells = outTime;
            obj.TempDataFileObj.ROImzVec = mzroi_end;
        end

        function obj = AlignScans(obj,modeFlag)
            mzQuan = obj.mzQuantil;
            mzEstim = obj.mzEstimMethod;
            mzCorr = obj.mzCorrectionMethod;
            switch modeFlag
                case "preview"
                    PeakCells = obj.TempDataFileObj.ROICells(1,1);
                otherwise
                    PeakCells = obj.TempDataFileObj.ROICells;
            end
            for id = 1:size(PeakCells,1)
                % perform Spectral Alignment
                [~, PeakCells{id,1}]= mspalign(PeakCells{id,1},'Quantile',mzQuan,'EstimationMethod',mzEstim,'CorrectionMethod',mzCorr,'ShowEstimation',false);
            end
            obj.TempDataFileObj.ROICells = PeakCells;
        end

        function obj = CorrectBaseline(obj,modeFlag)

            switch modeFlag
                case "preview"
                    MSroi = obj.TempDataFileObj.ROICells(1,1);
                    time=obj.TempDataFileObj.TimeCells(1,1);

                case "batch"
                    MSroi = obj.TempDataFileObj.ROICells;
                    time=obj.TempDataFileObj.TimeCells;
            end
            WSize = obj.WindowSize;
            SSize = obj.StepSize;
            RegMethod = obj.RegressionMethod;
            EstMethod =obj.EstimationMethod;
            SmooMethod = obj.SmoothMethod;
            Quan = obj.QuantilVal;
            parfor id = 1:size(MSroi,1)
                oldSize=size(MSroi{id,1});
                %depad Array
                MSroiTemp = MSroi{id,1};
                [MSroiTemp,timeTemp] = depadArrays(MSroiTemp,time{id,1});
                MSroiTemp = msbackadj(timeTemp,MSroiTemp,'WindowSize',WSize,'StepSize',SSize,'RegressionMethod',RegMethod,'EstimationMethod',EstMethod,'SmoothMethod',SmooMethod,'QuantileValue',Quan,'PreserveHeights',true);
                %remove negative, NaN and re-pad Array
                MSroiTemp=max(MSroiTemp,0);
                MSroiTemp(isnan(MSroiTemp))=0;
                [MSroi{id,1},time{id,1}] = repadArrays(MSroiTemp,timeTemp,oldSize);
                % set possible negative values to 0
                MSroi{id,1} = max(MSroi{id,1},0);
            end
            obj.TempDataFileObj.ROICells = MSroi;
            obj.TempDataFileObj.TimeCells = time;
        end

        function obj = SmoothPeaks(obj,modeFlag)
            switch modeFlag
                case "preview"
                    MSroi = obj.TempDataFileObj.ROICells(1,1);
                    time=obj.TempDataFileObj.TimeCells(1,1);
                case "batch"
                    MSroi = obj.TempDataFileObj.ROICells;
                    time=obj.TempDataFileObj.TimeCells;

            end
            Frame = obj.FrameSize;
            Deg = obj.Degree;
            parfor id = 1:size(MSroi,1)
                %depad Array
                oldSize=size(MSroi{id,1});
                MSroiTemp = MSroi{id,1};
                [MSroiTemp,timeTemp] = depadArrays(MSroiTemp,time{id,1});
                MSroiTemp = mssgolay(timeTemp,MSroiTemp,'Span',Frame,'Degree',Deg);
                %remove negative, NaN and re-pad Array
                MSroiTemp=max(MSroiTemp,0);
                MSroiTemp(isnan(MSroiTemp))=0;
                [MSroi{id,1},time{id,1}] = repadArrays(MSroiTemp,timeTemp,oldSize);
                % set possible negative values to 0
                MSroi{id,1} = max(MSroi{id,1},0);
            end
            obj.TempDataFileObj.ROICells = MSroi;
            obj.TempDataFileObj.TimeCells = time;
        end

        function obj = AlignPeaks(obj,modeFlag)

            switch modeFlag
                case "preview"
                    MSroi = obj.TempDataFileObj.ROICells(1,1);
                    time=obj.TempDataFileObj.TimeCells(1,1);
                case "batch"
                    MSroi = obj.TempDataFileObj.ROICells;
                    time=obj.TempDataFileObj.TimeCells;
            end
            WSR = obj.WindowSizeRatio;
            I = obj.Iterations;
            GS = obj.GridSteps;
            SS = obj.SearchSpace;
            ShiftVal = [obj.maxshiftneg*-1,obj.maxshiftpos];
            PW = obj.PulseWidth;
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

            parfor n=1:size(MSroi,2)
                ROI=reshape(MSroi(:,n),maxScan,[]);
                ROI=msalign(TimeVec,ROI,refTime(n),'MaxShift',ShiftVal,...
                    'WidthOfPulses',PW,'WindowSizeRatio',WSR,'Iterations',...
                    I,'GridSteps',GS,'SearchSpace',SS);
                ROI(isnan(ROI))=0; %remove possible NaN
                MStemp{1,n} = reshape(ROI,[],1);
            end
            MSroi = cell2mat(MStemp);
            obj.TempDataFileObj.ROICells = mat2cell(MSroi,splitVar);
            obj.TempDataFileObj.TimeCells = mat2cell(time,splitVar);
        end

        function obj = FinalizeROI(obj)
            time = obj.TempDataFileObj.TimeCells;
            MSroi = obj.TempDataFileObj.ROICells;
            PaddedSize = zeros(size(MSroi));
            timeTemp=horzcat(time{:});
            timeTemp(any(timeTemp==0,2),:)=[];
            obj.ScanFrequency=mean(diff(timeTemp),'all');
            MaxPW=round(obj.maxWidth*1.5/obj.ScanFrequency);
            parfor i = 1:size(MSroi,1)
                MStemp=MSroi{i,1};
                MStemp=max(MStemp,0);
                MStemp(isnan(MStemp))=0;
                MStemp=padarray(MStemp,MaxPW,0,'post');
                MSroi{i,1}=MStemp;
                time{i,1}=padarray(time{i,1},MaxPW,0,'post');
                PaddedSize(i) = length(time{i,1});
            end
            obj.nScansPadded = PaddedSize;
            obj.nScans(length(PaddedSize)+1:end) = [];
            MStemp = vertcat(MSroi{:});
            %remove empty columns
            id = all(MStemp < obj.thresh,1);
            MStemp(:,id) = [];
            MStemp=max(MStemp,0);
            if obj.BLKSubtraction == true
                obj.TempDataFileObj.ROIMatBLK(:,id) = [];
            end
            obj.TempDataFileObj.ROImzVec(id) = [];
            obj.TempDataFileObj.ROIMat = sparse(MStemp);
            obj.TempDataFileObj.timeVec = round(vertcat(time{:}),1);
        end

        function [PeakData,TimeData,PrecursorData,ColType,ColEnergy]= MS2CleanUp(obj,PeakData,TimeData,PrecursorData,ColType,ColEnergy)
            %remove empty scans and rescale intensities
            %% Clean Data
            for k = 1 : size(PeakData,1)
                Peak = PeakData{k,1};
                idx = cellfun(@isempty,Peak);
                Peak(idx,:) = [];
                TimeData{k,1}(idx,:) = [];
                PrecursorData{k,1}(idx,:) = [];
                ColType{k,1}(idx,:) = [];
                ColEnergy{k,1}(idx,:) = [];
                parfor n = 1:numel(Peak)
                    Peak{n,1}(:,2) = Peak{n,1}(:,2)/max(Peak{n,1}(:,2));
                end
                PeakData{k,1}=Peak;

            end
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