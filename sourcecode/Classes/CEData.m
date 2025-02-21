classdef CEData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        GroupName (1,1) string
        SeparationType (1,1) string = "CE"
    end

    methods
        function obj = CEData(groupNumber)
            if nargin == 0
                groupNumber = 0;
            end
            obj = obj@RawData;
            obj.GroupName = "Group " + groupNumber;
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
        
        %% helper functions
        
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

        function IntResults = CWTIntegrate(obj,Index)
            Mat = obj.TempDataFileObj.ROIMat;
            Mat = Mat(:,Index);
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
    
end