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
        scanPolarities      (1,1) string {mustBeMember(scanPolarities,["positive","negative","both"])} = "positive"
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
                    case "CDF"
                        [FileInfo(n),RetentionTimes{n},TIC{n},BPC{n},polarityCells{n}] = mzCDFinfo(FileLoc{n});
                end
            end
            obj.RawDataFileObj.PreviewTICs = TIC;
            obj.RawDataFileObj.PreviewBPCs = BPC;
            obj.RawDataFileObj.PreviewTimes = RetentionTimes;

            polarity = vertcat(polarityCells{:});

            test = strcmp(polarity,"+");
            if all(test)
                polarity = "positive";
            elseif all(~test)
                polarity = "negative";
            else
                polarity = "both";
            end
            obj.scanPolarities = polarity;

            obj.RawDataFileObj.DataMS1 ={[]};
            obj.RawDataFileObj.DataMS2 ={[]};

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

        function obj = ReadData(obj,DataLoc,separationType)
            nFiles = size(DataLoc,1);
            %preallocation
            ms1 = cell(nFiles,1);
            ms2 = cell(nFiles,1);
            parfor n = 1:nFiles
                %filetype check
                FileType = strsplit(DataLoc(n),'.');
                FileType = FileType(end);
                switch FileType
                    case "mzML"
                        [ms1{n,1},ms2{n,1}] = readmzML_MSandMS2(DataLoc{n});
                    case "mzXML"
                        [ms1{n,1},ms2{n,1}] = readmzXML_MSandMS2(DataLoc{n});
                    case "CDF"
                        ms1{n,1} =  readmzCDF(DataLoc{n});
                end
            end

            for n = 1:nFiles
                %remove possible empty scans
                emptyScans = cellfun(@isempty, ms1{n,1}.profileDataMS1);
                ms1{n,1}.profileDataMS1(emptyScans,:) = [];
                ms1{n,1}.timeDataMS1(emptyScans,:) = [];
                ms1{n,1}.polarityMS1(emptyScans,:) = [];

                emptyScans = cellfun(@isempty, ms2{n,1}.profileDataMS2);
                ms2{n,1}.profileDataMS2(emptyScans,:) = [];
                ms2{n,1}.timeDataMS2(emptyScans,:) = [];
                ms2{n,1}.polarityMS2(emptyScans,:) = [];
                ms2{n,1}.precursorMass(emptyScans,:) = [];
                ms2{n,1}.fragmentationEnergy(emptyScans,:) = [];
                ms2{n,1}.fragmentationType(emptyScans,:) = [];


                switch separationType
                    case "GC"
                        %centroid profile data
                        ms1{n,1}.centroidDataMS1 = CentroidScans(ms1{n,1}.profileDataMS1);

                    otherwise
                        %convert MS1 and MS2 precursor to molecular mass
                        ms1{n,1}.profileDataMS1 = ConvertScans2MolecularMass(ms1{n,1}.profileDataMS1,ms1{n,1}.polarityMS1);

                        precursors = ms2{n,1}.precursorMass;
                        modifier = ones(size(precursors))*1.007825;
                        idx = ms2{n,1}.polarityMS2 == "+";
                        modifier(idx) = modifier(idx)*-1;
                        ms2{n,1}.precursorMassCorrected = precursors + modifier;

                        %centroid profile data
                        ms1{n,1}.centroidDataMS1 = CentroidScans(ms1{n,1}.profileDataMS1);

                        %compress profile MS1 data
                        ms1{n,1} = DataCleanUp(ms1{n,1});

                        %compress MS2 data and store
                        ms2{n,1}.centroidDataMS2 = CentroidScans(ms2{n,1}.profileDataMS2);
                end
            end
            %store data
            obj.RawDataFileObj.DataMS1 = ms1;
            obj.RawDataFileObj.DataMS2 = ms2;

        end

        %% Data Processing

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

        function IntegrationData = AssignRT2SampleFile(obj,IntegrationData)
            test=cumsum(obj.nScansPadded)';
            nFiles = length(obj.nScansPadded);
            parfor f = 1:length(IntegrationData)
                fileID = zeros(size(IntegrationData(f).peakLocation));
                peakLocation = IntegrationData(f).peakLocation;
                for n = 1:numel(fileID)
                    val = peakLocation(n,1);
                    val = val < test;
                    val = sum(val,2)-1;
                    fileID(n,1) = abs(val-nFiles);
                end
                IntegrationData(f).fileID = fileID;
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

        function Output = GroupAndSampleScaling(obj,Output)
            if ~ Output.dataSize == 0
                features = Output.feature;
                for n = 1:length(features)
                    %GroupScale
                    features(n).peakHeights = features(n).peakHeights/obj.GroupScale;
                    features(n).peakAreas = features(n).peakAreas/obj.GroupScale;
                    %SampleScale
                    features(n).peakHeights = features(n).peakHeights./obj.SampScale;
                    features(n).peakAreas = features(n).peakAreas./obj.SampScale;
                end
                Output.feature = features;
            end
        end

        function obj = CutScansToSize(obj)
            StartTime = obj.Start;
            EndTime = obj.End;
            data = obj.RawDataFileObj.DataMS1;
            tempPeakData = cell(size(data));
            tempTimeData = tempPeakData;
            parfor n = 1:height(data)
                tempPeakData{n,1} = data{n,1}.centroidDataMS1;
                tempTimeData{n,1} = data{n,1}.timeDataMS1;
            end
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




        function IntResults = FilterPeaks(obj,IntResults,Noise)
            % Filters identified peaks from AutoCWT
            MinPWDataPoints=floor(obj.minWidth/obj.ScanFrequency);
            MaxPWDataPoints=ceil(obj.maxWidth/obj.ScanFrequency);
            maxSN = obj.minSignalNoise;


            parfor n = 1:length(IntResults)
                %check empty input
                if isempty(IntResults(n).peakLocation)
                    continue
                end

                %% Peak filter
                %remove duplicate peaks
                out = unique([IntResults(n).peakLocation,IntResults(n).peakStartLocation,IntResults(n).peakEndLocation,IntResults(n).peakHeight,],'rows','stable');
                IntResults(n).peakLocation = out(:,1);
                IntResults(n).peakStartLocation = out(:,2);
                IntResults(n).peakEndLocation = out(:,3);
                IntResults(n).peakHeight = out(:,4);

                %preallocate indexarray
                idx = false(size(IntResults(n).peakLocation));

                %remove peaks with wrong boundaries
                id = IntResults(n).peakStartLocation>=IntResults(n).peakEndLocation;
                idx = idx | id;

                %remove peaks with height = 0
                id = IntResults(n).peakHeight == 0;
                idx = idx | id;

                %remove peaks with bad Peak asymmetry
                symmetry = (IntResults(n).peakEndLocation - IntResults(n).peakLocation)./(IntResults(n).peakLocation - IntResults(n).peakStartLocation);
                id = symmetry<0.3 | symmetry>3;
                idx = idx | id;

                %less than minimum peak width
                id = IntResults(n).peakEndLocation-IntResults(n).peakStartLocation < MinPWDataPoints;
                IntResults(n).minWidthFiltered=sum(id);
                idx = idx | id;

                %more than maximum peak width
                id=IntResults(n).peakEndLocation - IntResults(n).peakStartLocation > MaxPWDataPoints;
                IntResults(n).maxWidthFiltered=sum(id);
                idx = idx | id;

                %S/N peak rejection
                IntResults(n).signal2Noise = IntResults(n).peakHeight ./ Noise(n);
                id = IntResults(n).signal2Noise < maxSN;
                IntResults(n).signal2NoiseFiltered = sum(id);
                idx = idx | id;

                % remove identified peaks
                IntResults(n).peakLocation(idx) = [];
                IntResults(n).peakStartLocation(idx) = [];
                IntResults(n).peakEndLocation(idx) = [];
                IntResults(n).peakHeight(idx) = [];
                IntResults(n).signal2Noise(idx) = [];
            end

            %calculate peak entropy and filter after first filter round
            % possible wrong peak bounderies causes errors
            IntResults = obj.CalculatePeakEntropy(IntResults);

            %determine entropy bins
            if obj.entropyFilter == true
                allEntropy = vertcat(IntResults(:).entropy);
                [~,binedges] = histcounts(allEntropy,'BinMethod','auto');
                switch obj.entropyStrength
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

            parfor n = 1:length(IntResults)
                %check empty input
                if isempty(IntResults(n).peakLocation)
                    continue
                end
                %entropy peak rejection
                id = IntResults(n).entropy > medianEntropy;
                IntResults(n).entropyFiltered = sum(id);

                IntResults(n).peakLocation(id) = [];
                IntResults(n).peakStartLocation(id) = [];
                IntResults(n).peakEndLocation(id) = [];
                IntResults(n).peakHeight(id) = [];
                IntResults(n).signal2Noise(id) = [];
            end

            %remove features without peaks
            id = false(length(IntResults),1);
            for n = 1:length(IntResults)
                id(n) = isempty(IntResults(n).peakLocation);
            end
            IntResults(id) = [];

        end


        function outputStruct = FindOriginalScans(obj,inputStruct)
            outputStruct = inputStruct;
            error = obj.mzerror;
            errorUnit = obj.mzErrorUnit;

            allScans = obj.RawDataFileObj.profileDataMS1;
            % append all scans with spacers in between, to match processing
            % indices
            for n = 1:numel(allScans)
                temp = allScans{n,1};
                temp(obj.nScansPadded(n),1) = {[]};
                allScans{n,1} = temp;
            end
            allScans = vertcat(allScans{:});
            numFiles = numel(obj.Files);

            parfor n = 1: length(inputStruct)
                avgSpectra = cell(1,numFiles);
                borders = inputStruct(n).peakBorders;

                for f = 1:numFiles
                    %check if borders contain NaN then skip iteration
                    if any(isnan(borders(:,f)))
                        continue
                    end
                    %select spectra in peak range
                    scans = allScans(borders(1,f):borders(2,f));
                    %remove possible empty scans
                    scans(cellfun(@isempty, scans)) = [];
                    if numel(scans) > 1 %average scan if multiple are present
                        times = 1:numel(scans);
                        [mzroi,MSroi,~] = ROIpeaks3(scans,0,error,errorUnit,1,times);
                        %calculate average spectrum
                        MSroi = mean(MSroi);
                        %rescale
                        MSroi = MSroi./max(MSroi,[],"all");
                        %reorder output
                        avgSpectra{1,f} = [mzroi;MSroi]';
                    else
                        avgSpectra{1,f} = scans{1,1};
                    end
                end
                outputStruct(n).spectrumMS1 = avgSpectra;
            end
        end


    end

    methods (Static)

        function IntegrationStruct = CalculatePeakEntropy(IntegrationStruct)
            % Calculates Peak entropy for all peaks
            for f = 1:length(IntegrationStruct)
                %check for no peaks, then skip iteration
                if isempty(IntegrationStruct(f).peakStartLocation)
                    continue
                end
                D = diff(IntegrationStruct(f).XIC(:,2));
                p = zeros(size(IntegrationStruct(f).peakLocation));
                for n = 1:numel(p)
                    %extract peak range
                    Peak = D(IntegrationStruct(f).peakStartLocation(n,:):IntegrationStruct(f).peakEndLocation(n,:));
                    maxidx = IntegrationStruct(f).peakLocation(n)-IntegrationStruct(f).peakStartLocation(n);
                    % check normal or variant point , variant point = 1
                    premax = Peak(1:maxidx-1)<0;
                    postmax = Peak(maxidx+1:end)>0;
                    VarPoints = [premax; false; postmax];
                    %calculate probability of variant point
                    p(n,1) = sum(VarPoints)/numel(VarPoints);
                end
                %calculate entropy
                PeakEntropy = -p.*log2(p)-(1-p).*log2(1-p);
                PeakEntropy(isnan(PeakEntropy)) = 0;
                %store values
                IntegrationStruct(f).entropy = PeakEntropy;
            end
        end

        function OutArray = FileSortPeaks(InArray)
            %% sort struct contents to original file

            %check if GC or LC/CE
            isGC = isscalar(InArray);

            %preallocate output Feature struct
            OutArray = InArray;

            nFiles = max(vertcat(InArray(:).fileID));
            %sort
            parfor n = 1:length(InArray)
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

                id = InArray(n).fileID;
                for fileID = 1:nFiles
                    % sort peaks into respective file cells
                    idx = id == fileID;
                    peakLcell{fileID} = InArray(n).peakLocation(idx);
                    peakRTcell{fileID} = InArray(n).peakRetentionTime(idx);
                    peakSLcell{fileID} = InArray(n).peakStartLocation(idx);
                    peakELcell{fileID} = InArray(n).peakEndLocation(idx);
                    peakHcell{fileID} = InArray(n).peakHeight(idx);
                    peakAcell{fileID} = InArray(n).peakArea(idx);
                    entroCell{fileID} = InArray(n).entropy(idx);
                    s2ncell{fileID} = InArray(n).signal2Noise(idx);
                    fileidcell{fileID} = InArray(n).fileID(idx);
                    if isGC
                        massCell{fileID} = InArray(n).mass(idx);
                        sMS2cell{fileID} = InArray(n).spectrumMS2(idx);
                    end
                end
                %store temp cells into output
                OutArray(n).peakLocation = peakLcell;
                OutArray(n).peakRetentionTime = peakRTcell;
                OutArray(n).peakStartLocation = peakSLcell;
                OutArray(n).peakEndLocation = peakELcell;
                OutArray(n).peakHeight = peakHcell;
                OutArray(n).peakArea = peakAcell;
                OutArray(n).entropy = entroCell;
                OutArray(n).signal2Noise = s2ncell;
                OutArray(n).fileID = fileidcell;

                if isGC
                    OutArray(n).mass = massCell;
                    OutArray(n).spectrumMS2 = sMS2cell;
                end
            end
        end


        function [PeakData,TimeData,PrecursorData,ColType,ColEnergy]= MS2CleanUp(PeakData,TimeData,PrecursorData,ColType,ColEnergy)
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

        function integrationStruct = FinalizeIntegrationOutput(integrationStruct,times)
            % Performs Integration of found Peaks and gathers retention times

            %check for empty struct
            if isempty(integrationStruct)
                return
            end

            parfor nfeats = 1:length(integrationStruct)
                EIC = full(integrationStruct(nfeats).XIC);
                areas = zeros(size(integrationStruct(nfeats).peakLocation));

                for n=1:numel(areas)
                    areas(n,1)=trapz(EIC(integrationStruct(nfeats).peakStartLocation(n,1):integrationStruct(nfeats).peakEndLocation(n,1)));
                end

                integrationStruct(nfeats).peakArea = areas;
                integrationStruct(nfeats).peakRetentionTime = times(integrationStruct(nfeats).peakLocation);
            end
        end


    end
end