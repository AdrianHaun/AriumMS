classdef GCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        GroupName (1,1) string
        SeparationType (1,1) string = "GC"
        EISpectra (:,1) cell
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
            obj.MinWidthFiltered = IntegrationData.minWidthFiltered;
            obj.MaxWidthFiltered = IntegrationData.maxWidthFiltered;
            obj.SNFiltered = IntegrationData.signal2NoiseFiltered;
            % calculate median entropy
            obj.MedianEntropy=median(IntegrationData.entropy,'omitnan');
            % Apply Entropy filter
            if obj.entropyFilter == true
                IntegrationData = obj.FilterbyEntropy(IntegrationData,obj.MedianEntropy);
                obj.EntropyFiltered = IntegrationData.entropyFiltered;
            end

            IntegrationData = obj.AssignRT2SampleFile(IntegrationData);

            % Build Storage Arrays and filter by number of occurences
            [Output,obj] = obj.BuildStorageArrays(IntegrationData);

            %check for empty Output
            if isempty(Output.dataSize)
                Output.dataSize = 0;
            end

            Output = obj.GroupAndSampleScaling(Output);

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
            %output preallocation
            IntResults = struct( ...
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
            %tics = mat2cell(tics,obj.nScansPadded);
            times = obj.TempDataFileObj.timeVec;
            IntResults.XIC = [tics,times];
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
            IntResults.peakLocation = peaks(:,1);
            IntResults.peakStartLocation = peaks(:,2);
            IntResults.peakEndLocation = peaks(:,3);
            IntResults.peakHeight = peaks(:,4);
            peaks = [];
            IntResults = obj.FilterPeaks(IntResults,MinPWDataPoints,MaxPWDataPoints,minSN,noise);
            % entropy calculation
            IntResults.entropy = CalculatePeakEntropy(IntResults,full(currentTIC));
            IntResults = FinalizeIntegrationOutput(IntResults,currentTIC,currentTime);
            IntResults = obj.GatherEISpectra(IntResults);
        end

        function [IntResults,obj] = GatherEISpectra(obj,IntResults)
            ROI = obj.TempDataFileObj.ROIMat;
            ROImz =  obj.TempDataFileObj.ROImzVec;
            FoundSpectra = cell(size(IntResults.peakLocation));
            MolecularMass = zeros(size(FoundSpectra));
            peakWidths = [IntResults.peakStartLocation,IntResults.peakEndLocation];
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
            IntResults.spectrumMS2 = FoundSpectra;
            IntResults.mass = MolecularMass;
        end

        function [output,obj] = BuildStorageArrays(obj,IntegrationResults,varargin)
            nFiles = numel(obj.Files);

            %preallocate Output struct
            output = struct(...
                "feature",[],...
                "minWidthFiltered",[],...
                "maxWidthFiltered",[],...
                "entropyFiltered",[],...
                "signal2NoiseFiltered",[],...
                "occurenceFiltered",[],...
                "groupName",string,...
                "fileNames",string,...
                "dataSize",[],...
                "separationType",string);

            %store group infos
            output.minWidthFiltered = IntegrationResults.minWidthFiltered;
            output.maxWidthFiltered = IntegrationResults.maxWidthFiltered;
            output.entropyFiltered = IntegrationResults.entropyFiltered;
            output.signal2NoiseFiltered = IntegrationResults.signal2NoiseFiltered;
            output.groupName = obj.GroupName;
            output.fileNames = obj.FileNames;
            output.separationType = obj.SeparationType;


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


            IntegrationResults = FileSortPeaks(IntegrationResults);

            %gather tolerances
            TimeTolerance = obj.RTTol;

            if isscalar(varargin)
                minDataPoints = nFiles;
                isISIntegration = true;
            else
                minDataPoints = ceil(nFiles*obj.minOccurence);
                isISIntegration = false;
            end

            %match features and store in feature struct
            uniqueFeatures = [vertcat(IntegrationResults.mass{:}),vertcat(IntegrationResults.peakRetentionTime{:})];
            uniqueFeatures = unique(uniqueFeatures,"rows");

            numFeatures = height(uniqueFeatures);

            %preallocat feature Storage
            featureStruct = repmat(featureStruct,numFeatures,1);


            parfor n = 1:numFeatures
                currentFeature = uniqueFeatures(n,:);
                % preallocate temp storages
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
                for file = 1:nFiles
                    idm = IntegrationResults.mass{file} == currentFeature(1,1);
                    idRT = abs(IntegrationResults.peakRetentionTime{file} - currentFeature(1,2)) <= TimeTolerance;
                    idx = idm & idRT;
                    if sum(idx) == 0 %no matching peaks
                        continue
                    elseif sum(idx) > 1 %split peak, ignore peak with lower intensity
                        h = max(IntegrationResults.peakHeight{file}(idx));
                        idx = idx & IntegrationResults.peakHeight{file} == h;
                    end
                    areas(1,file) = IntegrationResults.peakArea{file}(idx);
                    heights(1,file) = IntegrationResults.peakHeight{file}(idx);
                    retentionTimes(1,file) = IntegrationResults.peakRetentionTime{file}(idx);
                    peakLocation(1,file) = IntegrationResults.peakLocation{file}(idx);
                    peakBorders(1,file) = IntegrationResults.peakStartLocation{file}(idx);
                    peakBorders(2,file) = IntegrationResults.peakEndLocation{file}(idx);
                    signal2Noise(1,file) = IntegrationResults.signal2Noise{file}(idx);
                    entropy(1,file) = IntegrationResults.entropy{file}(idx);
                    spectrum{1,file} = IntegrationResults.spectrumMS2{file}(idx);
                    xic{1,file} = IntegrationResults.XIC(peakBorders(1,file):peakBorders(1,file),:);
                end
                %store matching features
                featureStruct(n).mass_measured = currentFeature(1,1);
                featureStruct(n).retentionTime = currentFeature(1,2);
                featureStruct(n).peakHeights = heights;
                featureStruct(n).peakAreas = areas;
                featureStruct(n).peakLocations = peakLocation;
                featureStruct(n).peakBorders = peakBorders;
                featureStruct(n).retentionTimes = retentionTimes;
                featureStruct(n).signal2Noise = signal2Noise;
                featureStruct(n).entropy = entropy;
                featureStruct(n).spectrumMS2 = spectrum;
                featureStruct(n).XIC = xic;
            end

            % remove features with less peaks than required minimum
            numElements = zeros(length(featureStruct),1);
            for ix = 1:length(featureStruct)
                numElements(ix) = nnz(~isnan(featureStruct(ix).peakHeights));
            end
            idx = numElements < minDataPoints;

            %sum number of removed peaks
            output.occurenceFiltered = sum(numElements(idx),"all");
            featureStruct(idx) = [];

            if isISIntegration == false
                %build featureID
                for ix = 1:length(featureStruct)
                    featureStruct(ix).featID = featureStruct(ix).mass_measured + "Da@" + featureStruct(ix).retentionTime + "s_" + obj.GroupName;
                end

                %finalize EI spectrum
                featureStruct = obj.FinalizeEISpectra(featureStruct);
                %gather original scans
                featureStruct = obj.FindOriginalScans(featureStruct);
            end

            output.feature = featureStruct;
            output.dataSize = length(output.feature);
            %local function
            function OutArray = FileSortPeaks(InArray)
                %% sort struct contents to original file and remove duplicate Features within one measurement

                %preallocate output Feature struct
                OutArray = struct( ...
                    "mass",cell(1), ...
                    "peakLocation",cell(1), ...
                    "peakRetentionTime",cell(1), ...
                    "peakStartLocation",cell(1), ...
                    "peakEndLocation",cell(1), ...
                    "peakHeight",cell(1), ...
                    "peakArea",cell(1), ...
                    "entropy",cell(1), ...
                    "signal2Noise",cell(1), ...
                    "spectrumMS2",cell(1), ...
                    "XIC",InArray.XIC, ...
                    "fileID",cell(1));

                %sort
                for fileID = 1:max(InArray.fileID)
                    id = InArray.fileID == fileID;
                    OutArray.mass{fileID} = InArray.mass(id);
                    OutArray.peakLocation{fileID} = InArray.peakLocation(id);
                    OutArray.peakRetentionTime{fileID} = InArray.peakRetentionTime(id);
                    OutArray.peakStartLocation{fileID} = InArray.peakStartLocation(id);
                    OutArray.peakEndLocation{fileID} = InArray.peakEndLocation(id);
                    OutArray.peakHeight{fileID} = InArray.peakHeight(id);
                    OutArray.peakArea{fileID} = InArray.peakArea(id);
                    OutArray.entropy{fileID} = InArray.entropy(id);
                    OutArray.signal2Noise{fileID} = InArray.signal2Noise(id);
                    OutArray.spectrumMS2{fileID} = InArray.spectrumMS2(id);
                    OutArray.fileID{fileID} = InArray.fileID(id);
                end
            end
        end

        function outputStruct = FinalizeEISpectra(obj,inputStruct)

            outputStruct = inputStruct;
            error = obj.mzerror;
            errorUnit = obj.mzErrorUnit;

            for n = 1:length(inputStruct)
                spectraCells = inputStruct(n).spectrumMS2;
                spectraCells = horzcat(spectraCells{:});
                if isscalar(spectraCells)
                    outputStruct(n).spectrumMS2 = spectraCells{:};
                else
                    %use ROI to sort values
                    %synthetic timevector
                    times = 1:numel(spectraCells);
                    [mzroi,MSroi,~] = ROIpeaks3(spectraCells',0,error,errorUnit,1,times);
                    %calculate average spectrum
                    MSroi = mean(MSroi);
                    %rescale
                    MSroi = MSroi./max(MSroi,[],"all");
                    %reorder output
                    outputStruct(n).spectrumMS2 = [mzroi;MSroi]';
                end
            end
        end

        function IntResults = FilterPeaks(obj,IntResults,MinPWDataPoints,MaxPWDataPoints,maxSN,Noise)
            % Filters identified peaks from AutoCWT
            %check empty input
            if isempty(IntResults.peakLocation)
                return
            end

            %% Peak filter
            %remove duplicate peaks
            out = unique([IntResults.peakLocation,IntResults.peakStartLocation,IntResults.peakEndLocation,IntResults.peakHeight,],'rows','stable');
            IntResults.peakLocation = out(:,1);
            IntResults.peakStartLocation = out(:,2);
            IntResults.peakEndLocation = out(:,3);
            IntResults.peakHeight = out(:,4);

            %preallocate indexarray
            idx = false(size(IntResults.peakLocation));

            %remove peaks with wrong boundaries
            id = IntResults.peakStartLocation>=IntResults.peakEndLocation;
            idx = idx | id;

            %remove peaks with height = 0
            id = IntResults.peakHeight == 0;
            idx = idx | id;

            %remove peaks with bad Peak asymmetry

            symmetry = (IntResults.peakEndLocation - IntResults.peakLocation)./(IntResults.peakLocation - IntResults.peakStartLocation);
            id = symmetry<0.3 | symmetry>3;
            idx = idx | id;

            %less than minimum peak width
            id = IntResults.peakEndLocation-IntResults.peakStartLocation < MinPWDataPoints;
            IntResults.minWidthFiltered=sum(id);
            idx = idx | id;

            %more than maximum peak width
            id=IntResults.peakEndLocation - IntResults.peakStartLocation > MaxPWDataPoints;
            IntResults.maxWidthFiltered=sum(id);
            idx = idx | id;

            %S/N peak rejection
            IntResults.signal2Noise = IntResults.peakHeight ./ Noise;
            id = IntResults.signal2Noise < maxSN;
            IntResults.signal2NoiseFiltered = sum(id);
            idx = idx | id;

            % remove identified peaks
            IntResults.peakLocation(idx) = [];
            IntResults.peakStartLocation(idx) = [];
            IntResults.peakEndLocation(idx) = [];
            IntResults.peakHeight(idx) = [];
            IntResults.signal2Noise(idx) = [];
        end
    end
end