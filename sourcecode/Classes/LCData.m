classdef LCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        GroupName (1,1) string
        SeparationType (1,1) string = "LC"
    end

    methods
        function obj = LCData(groupNumber)
            %Construct an instance of this class
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
            for n=1:nFiles
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
                    case "CDF"
                        [peakTemp,timeTemp] = mzcdf2peaks(mzcdfread(DataLoc{n},'Verbose',false));
                end
                %remove possible empty scans
                emptyScans = cellfun(@isempty, peakTemp);
                peakTemp(emptyScans) = [];
                timeTemp(emptyScans) = [];
                polarities{n}(emptyScans) = [];
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

            %Common Contaminant filter
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
            IntegrationData = obj.LCIntegrate(IDX);

            IntegrationData = obj.AssignRT2SampleFile(IntegrationData);
            IntegrationData = obj.FileSortPeaks(IntegrationData);
            IntegrationData = obj.mergeDuplicatePeaksWithinFile_LC(IntegrationData);

            %%%%%%%
            % % remove adducts
            % if obj.AdductFilter == true
            %     [IntegrationData,obj] = obj.FilterAdducts(IntegrationData);
            % end
            %%%%%%%
            

            % Build Storage Arrays and filter by number of occurences
            [Output,obj] = obj.BuildStorageArrays_LC(IntegrationData);
            % apply scaling
            Output = obj.GroupAndSampleScaling(Output);

            %gather MS2 spectra
            Output = obj.GatherFragmentSpectra(Output);

            %confirm same feature by MS2 comparison

           
            obj.Output = Output;

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
            ISIntegrationData = obj.LCIntegrate(ISid);
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
            ISData = obj.BuildStorageArrays_LC(ISIntegrationData,obj.ISMassFound);

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

        function IntResults = LCIntegrate(obj,Index)
            %gather data            
            Mat = obj.TempDataFileObj.ROIMat;
            Mat = Mat(:,Index);
            mzVec = obj.TempDataFileObj.ROImzVec;
            times = obj.TempDataFileObj.timeVec;

            % calculate EIC derivatives and store as sparse
            smoothed = smoothdata(Mat,"gaussian","omitnan","SmoothingFactor",0.1);
            Diff2 = zeros(length(times),size(Mat,2));
            Diff2(1:end-2,:) = diff(smoothed,2);
            numEIC = size(Mat,2);

            % prepare wavelet filter-bank
            FilterBank = cwtfilterbank("SignalLength",size(Diff2,1), ...
                "WaveletParameters",[3 4], ...
                "VoicesPerOctave",8, ...
                "SamplingPeriod",seconds(obj.ScanFrequency), ...
                "PeriodLimits",[seconds(obj.minWidth) seconds(obj.maxWidth)]);
            
            %preallocate storage struct
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
            IntResults = repmat(IntResults,numEIC,1);

            parfor id = 1:numEIC
                peaks = AutoCWT(Diff2(:,id),smoothed(:,id),FilterBank);
                eic = Mat(:,id);
                % Correct Peak Borders
                peaks = CWTBorderCorrection(peaks,eic,smoothed(:,id));
                IntResults(id).mass = mzVec(id);
                IntResults(id).peakLocation = peaks(:,1);
                IntResults(id).peakStartLocation = peaks(:,2);
                IntResults(id).peakEndLocation = peaks(:,3);
                IntResults(id).peakHeight = peaks(:,4);
                %store EIC
                IntResults(id).XIC = eic;
            end
            %filtere found peaks
            noise = std(Mat-smoothed);
            IntResults = obj.FilterPeaks(IntResults,noise);
            IntResults = obj.FinalizeIntegrationOutput(IntResults,times);
        end

        function outArray = mergeDuplicatePeaksWithinFile_LC(obj,inArray)

            rtTol = obj.RTTol;
            outArray = inArray;

            for n = 1:length(inArray) %sample loop
                
                for file = 1:numel(inArray(n).peakRetentionTime)
                    

                end
            end
        end

        function [output,obj] = BuildStorageArrays_LC(obj,IntegrationResults,varargin)
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
            output.minWidthFiltered = IntegrationResults(1).minWidthFiltered;
            output.maxWidthFiltered = IntegrationResults(1).maxWidthFiltered;
            output.entropyFiltered = IntegrationResults(1).entropyFiltered;
            output.signal2NoiseFiltered = IntegrationResults(1).signal2NoiseFiltered;
            output.fileNames = obj.FileNames;

            XIC = IntegrationResults(1).XIC;

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
            TimeTolerance = obj.RTTol;
            mzTolerance = obj.mzTol;
            mztolUnit = obj.mzTolUnit;

            if isscalar(varargin)
                minDataPoints = nFiles;
                isISIntegration = true;
            else
                minDataPoints = ceil(nFiles*obj.minOccurence);
                isISIntegration = false;
            end

            %match features and store in feature struct

            uniqueFeatures = [vertcat(IntegrationResults(:).mass),vertcat(IntegrationResults(:).peakRetentionTime)];

            % over preallocat feature Storage
            featureStruct = repmat(featureStruct,height(uniqueFeatures),1);
            n = 0;
            while ~isempty(uniqueFeatures)
                n = n+1; disp(n)
                currentFeature = uniqueFeatures(1,:);
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
                    switch mztolUnit
                        case "Da"
                            idm = abs(IntegrationResults(file).mass-currentFeature(1,1)) <= mzTolerance;
                        case "ppm"
                            idm = abs(IntegrationResults(file).mass-currentFeature(1,1))./currentFeature(1,1)*10^6 <= mzTolerance;
                    end
                    idRT = abs(IntegrationResults(file).peakRetentionTime - currentFeature(1,2)) <= TimeTolerance;
                    idx = idm & idRT;
                    % handle matching peaks
                    if sum(idx) == 0 %no matching peaks
                        continue
                    elseif sum(idx) > 1 %use peak with lower time tolerance
                        [~,idmin] = min(abs(IntegrationResults(file).peakRetentionTime - currentFeature(1,2)));
                        idx = false(size(idx));
                        idx(idmin) = true;
                    end
                    %store found peak information
                    areas(1,file) = IntegrationResults(file).peakArea(idx);
                    heights(1,file) = IntegrationResults(file).peakHeight(idx);
                    retentionTimes(1,file) = IntegrationResults(file).peakRetentionTime(idx);
                    peakLocation(1,file) = IntegrationResults(file).peakLocation(idx);
                    peakBorders(1,file) = IntegrationResults(file).peakStartLocation(idx);
                    peakBorders(2,file) = IntegrationResults(file).peakEndLocation(idx);
                    signal2Noise(1,file) = IntegrationResults(file).signal2Noise(idx);
                    entropy(1,file) = IntegrationResults(file).entropy(idx);
                    spectrum{1,file} = IntegrationResults(file).spectrumMS2(idx);
                    xic{1,file} = XIC(peakBorders(1,file):peakBorders(1,file),:);

                    %delete peaks from input struct
                    IntegrationResults(file).mass(idx) = [];
                    IntegrationResults(file).peakArea(idx) = [];
                    IntegrationResults(file).peakHeight(idx) = [];
                    IntegrationResults(file).peakRetentionTime(idx) = [];
                    IntegrationResults(file).peakLocation(idx) = [];
                    IntegrationResults(file).peakStartLocation(idx) = [];
                    IntegrationResults(file).peakEndLocation(idx) = [];
                    IntegrationResults(file).signal2Noise(idx) = [];
                    IntegrationResults(file).entropy(idx) = [];
                    IntegrationResults(file).spectrumMS2(idx) = [];


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

                %update remaining features
                uniqueFeatures(1,:) = [];
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
                    featureStruct(ix).featID = featureStruct(ix).mass_measured + "Da@" + featureStruct(ix).retentionTime + "s_";
                end

                %gather original scans
                featureStruct = obj.FindOriginalScans(featureStruct);
            end

            output.feature = featureStruct;
            output.dataSize = length(output.feature);
        end

    end
end