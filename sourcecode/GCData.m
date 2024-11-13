classdef GCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        GroupName (1,1) string
        SeparationType (1,1) string = "GC"
    end

    methods
        function obj = GCData(groupNumber)
            %Construct an instance of this class
             if nargin == 0
                groupNumber = 0;
            end
            obj = obj@RawData;
            obj.GroupName = "Group " + groupNumber;
            % set default parameters
            obj.mzerror = 0.1;
            obj.mzErrorUnit = "Da";
            obj.minroi = 10;
            obj.minWidth = 0.8;
            obj.maxWidth = 10;
            obj.mzTol = 0.05;
            obj.mzTolUnit = "Da";
        end


        function obj = ReadData(obj,DataLoc,~)
            nFiles = size(DataLoc,1);
            %preallocation
            Peaks=cell(nFiles,1);
            times=cell(nFiles,1);
            fileType = obj.MSFileType;

            parfor n=1:nFiles
                peakTemp = [];
                timeTemp = [];
                %filetype check
                FileType=strsplit(DataLoc(n),'.');
                FileType=FileType(end);
                switch FileType
                    case "mzML"
                        [peakTemp,timeTemp] = readmzML(DataLoc{n},MSLevel=1);
                    case "mzXML"
                        [peakTemp,timeTemp] = readmzXML(DataLoc{n},MSLevel=1);
                end

                % when profile data then centroid scans
                if fileType == "profile"
                    peakTemp = CentroidScans(peakTemp);
                else
                    [peakTemp,timeTemp] = DataCleanUp(peakTemp,timeTemp);
                end
                Peaks{n,1} = peakTemp;
                times{n,1} = timeTemp;
            end
                obj.RawDataFileObj.TimeDataMS1 = times;
                obj.RawDataFileObj.PeakDataMS1 = Peaks;
        end

        %% Data Processing
        function [Output,obj] = BatchProcess(obj,varargin)
            %check if old results exist and delete them
            if isfile(obj.ROIDataFile)
                delete(obj.ROIDataFile)
            end
            % check for OptimizationMode
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
            IntegrationData = obj.IntegrateGC;
            %Calculate number of removed features
            obj.MinWidthFiltered = sum(vertcat(IntegrationData{5,:}),"all");
            obj.MaxWidthFiltered = sum(vertcat(IntegrationData{6,:}),"all");
            obj.SNFiltered = sum(vertcat(IntegrationData{7,:}),"all");
            IntegrationData(5:7,:) = [];
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

        %% helper functions
        function IntResults = IntegrateGC(obj)

            %prepare TIC Data
            tics = sum(obj.TempDataFileObj.ROIMat,2);
            %tics = mat2cell(tics,obj.nScansPadded);
            times = obj.TempDataFileObj.timeVec;
            %gather parameters
            minSN = obj.minSignalNoise;
            MinPWDataPoints=floor(obj.minWidth/obj.ScanFrequency);
            MaxPWDataPoints=ceil(obj.maxWidth/obj.ScanFrequency);
            currentTIC = full(tics);
            currentTime = full(times);

            %calculate noise level
            smoothedTIC = smoothdata(currentTIC,"gaussian",4,"omitnan");
            noise = mean(std(currentTIC-smoothedTIC));
            currentTIC = sum(currentTIC,2);
            smoothedTIC = sum(smoothedTIC,2);

            [~,peakLoc,peakWidth] = findpeaks(currentTIC,"WidthReference","halfheight");
            %calculate initial borders and bring in correct form
            lowerBorders = floor(peakLoc-peakWidth/2);
            upperBorders = ceil(peakLoc+peakWidth/2);
            peaks = [peakLoc,lowerBorders,upperBorders];
            peaks = CWTBorderCorrection(peaks,currentTIC,smoothedTIC);
            [peaks,tempStorage] = FilterPeaks(peaks,MinPWDataPoints,MaxPWDataPoints,minSN,noise,currentTIC);
            IntResults = FinalizeIntegrationOutput(peaks,tempStorage,currentTIC,currentTime);
            IntResults = obj.GatherEISpectra(IntResults);
        end

        function [IntResults,obj] = GatherEISpectra(obj,IntResults)
            ROI = obj.TempDataFileObj.ROIMat;
            ROImz =  obj.TempDataFileObj.ROImzVec;
            FoundSpectra = cell(size(IntResults{1,1}));
            MolecularMass = zeros(size(FoundSpectra));
            peakWidths = IntResults{2,1}(:,2:3);
            EIlosses = load("MassListData.mat","EICommonLoss");
            EIlosses = EIlosses.EICommonLoss;
            parfor n = 1:height(MolecularMass)
                peakBorders = peakWidths(n,:);
                Spectras = ROI(peakBorders(1):peakBorders(2),:);
                %mean spectra
                Spectras = mean(Spectras);
                %normalize Spectras
                Spectras = Spectras./max(Spectras,[],"all");
                Spectras = [ROImz;full(Spectras)]';
                Spectras(Spectras(:,2)<0.05,:) = [];
                FoundSpectra{n,1} = Spectras;
                %identify molecular mass
                masses = flip(Spectras(:,1));
                hasMolecularMass = false;
                counter = 0;
                while hasMolecularMass == false & counter < numel(masses)
                    counter = counter+1;
                    possibleFragment = masses(counter)-EIlosses;
                    hasMolecularMass = any(min(abs(masses-possibleFragment'))<0.3);
                end
                if hasMolecularMass
                    MolecularMass(n,1) = masses(counter,1);
                else
                    MolecularMass(n,1) = Spectras(end,1);
                end
            end
            IntResults{9,1} = FoundSpectra;
            IntResults{10,1} = MolecularMass;
        end

        function [Output,obj] = BuildStorageArrays(obj,IntegrationResults,varargin)
            TimeTolerance = obj.RTTol;
            nFiles = size(obj.nScans,1);
            if isscalar(varargin)
                mzVector = varargin{1};
                minDataPoints = nFiles;
                isISIntegration = true;
            else
                mzVector = IntegrationResults{7,1};
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
    end
end