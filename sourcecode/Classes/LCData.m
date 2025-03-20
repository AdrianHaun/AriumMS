classdef LCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        SeparationType (1,1) string = "LC"
    end

    methods
        function obj = LCData(groupName,window)
            %Construct an instance of this class
            if nargin == 0
                groupName = 0;
                window = 0;
            end
            obj = obj@RawData(groupName,window);
        end

        %% Data Processing
        function [Output,obj] = BatchProcess(obj,varargin)
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
            else
                title = "Processing " + obj.GroupName;
                progressBar = uiprogressdlg(obj.mainWindow,"Title",title,"Message","Preparation",Value=0);
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

            progressBar.Message = "Loading files";
            %check if files already loaded then skip loading stage
            test = obj.RawDataFileObj.centroidedDataMS1;
            if isempty(test{1,1}) || size([obj.Files;obj.BlankFiles],1) ~= height(test)
                obj = obj.ReadData(FileLocs,obj.SeparationType);
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
            progressBar.Value = 0.33;
            tempp = [obj.TempDataFileObj.ROICells,obj.TempDataFileObj.TimeCells];
            assignin("base","ScansPreIsotope",tempp)

            % remove isotopes
            if obj.IsotopeFilter == true
                obj = obj.FilterIsotopesScanStage;
            end

            tempp = [obj.TempDataFileObj.ROICells,obj.TempDataFileObj.TimeCells];
            assignin("base","ScansPostIsotope",tempp)

            tempp = [obj.RawDataFileObj.profileDataMS1,obj.RawDataFileObj.timeDataMS1];
            assignin("base","ScansProfile",tempp)

            obj.nScans = cellfun(@numel,obj.TempDataFileObj.TimeCells);
            if obj.MSalign == true
                progressBar.Message = "Aligning MS Scans";
                obj = obj.AlignScans("batch");
                progressBar.Value = 0.4;
            end

            % ROI Search
            progressBar.Message = "Searching for ROIs";
            obj = obj.AutoROI("batch");
            progressBar.Value = 0.5;

            % Average BLK
            if obj.BLKSubtraction == true && nBLK > 1
                obj = obj.AverageBLK(nBLK);
                nData = size(obj.TempDataFileObj.ROICells,1); % update number of matrices
            end

            %Common Contaminant filter
            if obj.ContaminantFilter == true
                progressBar.Message = "Removing Contaminants";
                obj = obj.removeContaminants;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Baseline Correction
            if obj.BaseCorr == true
                progressBar.Message = "Correcting Baseline";
                obj = obj.CorrectBaseline("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Smoothing
            if obj.Smoothing == true
                progressBar.Message = "Smoothing Peaks";
                obj = obj.SmoothPeaks("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Peak Align
            if obj.Peakalign == true && nData > 1
                progressBar.Message = "Aligning Peaks";
                obj = obj.AlignPeaks("batch");
                progressBar.Value = progressBar.Value + 0.05;
            end

            if obj.BLKSubtraction == true % Separate Blank data from Sample data
                tempBLK = obj.TempDataFileObj.ROICells(end,1);
                obj.TempDataFileObj.ROIMatBLK=sparse(tempBLK{:});
                obj.TempDataFileObj.ROICells(end)=[];
                obj.TempDataFileObj.TimeCells(end)=[];
            end

            % subtract blank before IS normalization
            if obj.BLKSubtraction == true && obj.ISOrder == "BlankIS"
                progressBar.Message = "Subtracting Blank";
                peakCells = obj.TempDataFileObj.ROICells;
                BLKMat = obj.TempDataFileObj.ROIMatBLK;
                parfor id=1:size(peakCells,1)
                    peakCells{id,1}=peakCells{id,1}-BLKMat;
                    % set possible negative values to 0
                    peakCells{id,1} = max(peakCells{id,1},0);
                end
                obj.TempDataFileObj.ROICells = peakCells;
                progressBar.Value = progressBar.Value + 0.05;
            end
            
            % pad arrays with Maximum peak width*3 Scans to eliminate
            % integration interference between matrices
            obj = obj.FinalizeROI;

            clearvars -except obj progressBar
            %% Integration Stage
            % Find and Integrate IS separate
            if obj.ISTDCorr == true
                progressBar.Message = "Searching for Internal Standard";
                obj = obj.IntegrateIS;
                if ~isempty(obj.ISValue)
                    obj = obj.ISNormalize;
                end
                 progressBar.Value = progressBar.Value + 0.05;
            end
            % BLK Subtraction after IS Correction
            if obj.BLKSubtraction == true && obj.ISOrder == "ISBlank"
                progressBar.Message = "Subtracting Blank";
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
                progressBar.Value = progressBar.Value + 0.05;
            end
            % mass correction
            if obj.MassCal == true && ~isempty(obj.ISValue)
                progressBar.Message = "Performing IS mass correction";
                obj = obj.ISMassCorrection;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Integrate all Peaks
            IDX = true(1,size(obj.TempDataFileObj.ROIMat,2));
            progressBar.Message = "Integrating Peaks";
            IntegrationData = obj.LCIntegrate(IDX);
            progressBar.Value = 0.9;

            progressBar.Message = "Processing found Features";
            IntegrationData = obj.AssignRT2SampleFile(IntegrationData);
            IntegrationData = obj.FileSortPeaks(IntegrationData);

            %%%%%%%
            % % remove adducts
            % if obj.AdductFilter == true
            %     [IntegrationData,obj] = obj.FilterAdducts(IntegrationData);
            % end
            %%%%%%%
            

            % Build Storage Arrays and filter by number of occurences
            [Output,obj] = obj.BuildStorageArrays_LC(IntegrationData);

            %
            % Output = obj.ConfirmSameFeatureByIsotopeDistribution(Output);

            % Occurence filter
            Output = obj.OccurenceFilterFeatures(Output);

            %gather MS2 spectra
            Output = obj.GatherMS2Spectra(Output);
            
            %fill remaining fields
            Output = obj.FinalizeBatchOutput(Output);

            % apply scaling
            progressBar.Message = "Apply scaling";
            Output = obj.GroupAndSampleScaling(Output);
            
            progressBar.Message = "Group processing successful";
            progressBar.Value = 1;

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
            close(progressBar)
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
                IntResults(id).XIC = [times,eic];
            end
            %filtere found peaks
            noise = std(Mat-smoothed);
            IntResults = obj.FilterPeaks(IntResults,noise);
            IntResults = obj.FinalizeIntegrationOutput(IntResults,times);
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
            output.minWidthFiltered = sum(vertcat(IntegrationResults(:).minWidthFiltered));
            output.maxWidthFiltered = sum(vertcat(IntegrationResults(:).maxWidthFiltered));
            output.entropyFiltered = sum(vertcat(IntegrationResults(:).entropyFiltered));
            output.signal2NoiseFiltered = sum(vertcat(IntegrationResults(:).signal2NoiseFiltered));
            output.fileNames = obj.FileNames;
            output.groupName = obj.GroupName;
            output.separationType = obj.SeparationType;

            %remove unnecessary fields from input struct
            IntegrationResults = rmfield(IntegrationResults,["minWidthFiltered","maxWidthFiltered","entropyFiltered","signal2NoiseFiltered"]);

            emptyStruct = struct(...
                "featID",strings,...
                "mass_measured",[],...
                "retentionTime",NaN,...
                "adductType",strings,...
                "mass_corrected",[],...
                "formula",strings,...
                "peakHeights",NaN(1,nFiles),...
                "peakAreas",NaN(1,nFiles),...
                "peakLocations",NaN(1,nFiles),...
                "peakBorders",NaN(2,nFiles),...
                "retentionTimes",NaN(1,nFiles),...
                "signal2Noise",NaN(1,nFiles),...
                "entropy",NaN(1,nFiles),...
                "XIC",cell(1),...
                "spectrumMS1",cell(1),...
                "spectrumMS2",cell(1),...
                "asymmetry",[]);

            storedFeatures = cell(length(IntegrationResults),1);

            %gather tolerances
            timeTolerance = obj.RTTol;

            if isscalar(varargin)
                isISIntegration = true;
            else
                isISIntegration = false;
            end

            %match features and store in feature struct

            parfor featMass = 1:length(IntegrationResults)
                currentFeatureStruct = emptyStruct;
                currentFeatureStruct.mass_measured = IntegrationResults(featMass).mass;
                currentFeatureStruct.XIC = IntegrationResults(featMass).XIC;

                numPeaks = numel(vertcat(IntegrationResults(featMass).peakLocation{:}));

                %unpack data
                peakData = zeros(numPeaks,10);
                peakData(:,1) = vertcat(IntegrationResults(featMass).peakLocation{:});
                peakData(:,2) = vertcat(IntegrationResults(featMass).peakRetentionTime{:});
                peakData(:,3) = vertcat(IntegrationResults(featMass).peakStartLocation{:});
                peakData(:,4) = vertcat(IntegrationResults(featMass).peakEndLocation{:});
                peakData(:,5) = vertcat(IntegrationResults(featMass).peakHeight{:});
                peakData(:,6) = vertcat(IntegrationResults(featMass).peakArea{:});
                peakData(:,7) = vertcat(IntegrationResults(featMass).entropy{:});
                peakData(:,8) = vertcat(IntegrationResults(featMass).signal2Noise{:});
                peakData(:,9) = vertcat(IntegrationResults(featMass).fileID{:});
                peakData(:,10) = (peakData(:,4)-peakData(:,1))./(peakData(:,1)-peakData(:,3)); %asymmetry factor
                
                % preallocat current feature Storage
                currentFeatureStruct = repmat(currentFeatureStruct,height(peakData),1);
                
                %store first new entry
                currentFile = peakData(1,9);
                currentFeatureStruct(1).peakLocations(currentFile) = peakData(1,1);
                currentFeatureStruct(1).retentionTimes(currentFile) = peakData(1,2);
                currentFeatureStruct(1).retentionTime = peakData(1,2);
                currentFeatureStruct(1).peakBorders(:,currentFile) = [peakData(1,3);peakData(1,4)];
                currentFeatureStruct(1).peakHeights(currentFile) = peakData(1,5);
                currentFeatureStruct(1).peakAreas(currentFile) = peakData(1,6);
                currentFeatureStruct(1).entropy(currentFile) = peakData(1,7);
                currentFeatureStruct(1).signal2Noise(currentFile) = peakData(1,8);
                currentFeatureStruct(1).asymmetry = peakData(1,10);
                peakData(1,:) = [];

                %assign remaining peaks to features
                while ~isempty(peakData) 

                    currentRT = peakData(1,2);
                    currentFile = peakData(1,9);
                    currentAsymmetry = peakData(1,10);
                    id = abs(vertcat(currentFeatureStruct(:).retentionTime)-currentRT)<=timeTolerance;
                    matchingRT = sum(id);
                    
                    if matchingRT == 0 %no matching RT -> new Feature
                        currentFeatureStruct = storeInNewFeat(currentFeatureStruct,currentFile,peakData);

                    elseif matchingRT == 1 %single feature -> store
                        if isnan(currentFeatureStruct(id).peakLocations(currentFile)) %check if a peak is already present
                            currentFeatureStruct(id).peakLocations(currentFile) = peakData(1,1);
                            currentFeatureStruct(id).retentionTimes(currentFile) = peakData(1,2);
                            currentFeatureStruct(id).peakBorders(:,currentFile) = [peakData(1,3);peakData(1,4)];
                            currentFeatureStruct(id).peakHeights(currentFile) = peakData(1,5);
                            currentFeatureStruct(id).peakAreas(currentFile) = peakData(1,6);
                            currentFeatureStruct(id).entropy(currentFile) = peakData(1,7);
                            currentFeatureStruct(id).signal2Noise(currentFile) = peakData(1,8);
                            %average retentionTime
                            currentFeatureStruct(id).retentionTime = mean([currentFeatureStruct(id).retentionTime;peakData(1,2)],'omitnan');
                        else
                            currentFeatureStruct = storeInNewFeat(currentFeatureStruct,currentFile,peakData);
                        end

                    else %multiple matching features -> store based on asymmetry factor
                        [~,idA] = min(vertcat(currentFeatureStruct(id).asymmetry)-currentAsymmetry);
                        if isnan(currentFeatureStruct(idA).peakLocations(currentFile)) %check if a peak is already present
                            currentFeatureStruct(idA).peakLocations(currentFile) = peakData(1,1);
                            currentFeatureStruct(idA).retentionTimes(currentFile) = peakData(1,2);
                            currentFeatureStruct(idA).peakBorders(:,currentFile) = [peakData(1,3);peakData(1,4)];
                            currentFeatureStruct(idA).peakHeights(currentFile) = peakData(1,5);
                            currentFeatureStruct(idA).peakAreas(currentFile) = peakData(1,6);
                            currentFeatureStruct(idA).entropy(currentFile) = peakData(1,7);
                            currentFeatureStruct(idA).signal2Noise(currentFile) = peakData(1,8);
                            %average retentionTime
                            currentFeatureStruct(idA).retentionTime = mean([currentFeatureStruct(idA).retentionTime;peakData(1,2)],'omitnan');
                        else
                            currentFeatureStruct = storeInNewFeat(currentFeatureStruct,currentFile,peakData);
                        end
                    end
                    %remove stored peak from list
                    peakData(1,:) = [];
                end

                %remove empty structs
                id = isnan([currentFeatureStruct(:).retentionTime])';
                currentFeatureStruct(id) = [];

                %store currentFeatureStruct
                storedFeatures{featMass,1} = currentFeatureStruct;
            end

            %unzip features
            storedFeatures = vertcat(storedFeatures{:});
            %remove asymmetry field
            storedFeatures = rmfield(storedFeatures,"asymmetry");
            if isISIntegration == false
                %gather original scans
                storedFeatures = obj.FindOriginalScans(storedFeatures);
            end

            output.feature = storedFeatures;
        end

        function outputStruct = GatherMS2Spectra(obj,outputStruct)
            %check if MSn data is already loaded
            if isscalar(obj.RawDataFileObj.centroidedDataMS2)
                obj = obj.ReadData(obj.Files,obj.SeparationType);
            end

            %%%%%
            % test tolerances
            % rttol = obj.RTTol;
            % mztol = obj.mzTol;
            % mztolUnit = obj.mzTolUnit;

            mzTol = 0.05;
            mztolUnit = "Da";
            rttol = 10;
            %%%%%%
            

            times = obj.RawDataFileObj.timeDataMS2;
            times = vertcat(times{:});
            scans = obj.RawDataFileObj.centroidedDataMS2;
            scans = vertcat(scans{:});
            precursor = obj.RawDataFileObj.molecularPrecursorMass;
            precursor = vertcat(precursor{:});
                
            features = outputStruct.feature;

            parfor n = 1:length(features)
                idM = [];
                switch mztolUnit
                    case "Da"
                        idM = abs(precursor-features(n).mass_measured) <= mzTol;
                    case "ppm"
                        idM = abs(precursor-features(n).mass_measured)./features(n).mass_measured*10^6 <= mzTol;
                end
                idT = abs(times-features(n).retentionTime) <= rttol;
                id = idT & idM;
                foundScans = scans(id);
                %remove possible empty scans
                    foundScans(cellfun(@isempty, foundScans)) = [];
                    if numel(foundScans) >= 1 
                        foundScans = AlignSpectra(foundScans,"average","low");
                    else % no found scan
                        foundScans = [];
                    end
                    features(n).spectrumMS2 = foundScans;
            end
            outputStruct.feature = features;
        end
    end
end