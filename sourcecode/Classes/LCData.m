classdef LCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        separationType (1,1) string = "LC"
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

            fileArray = obj.dataFile;
            nBlanks = 0;
            if obj.useBlankSubtraction == true
                nBlanks = size(obj.blankFile,1);
                fileArray=[fileArray;obj.blankFile];
            end
            %remove possible empty
            iFile = cellfun(@isempty,fileArray);
            fileArray(iFile) = [];
            nData = numel(fileArray);

            progressBar.Message = "Loading files";
            %check if files already loaded then skip loading stage
            test = obj.RawDataFileObj.centroidDataMS1;
            if isempty(test{1,1}) || size([obj.fileName;obj.blankFile],1) ~= height(test)
                obj = obj.readData(fileArray,obj.separationType);
            end

            clearvars test fileArray id

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
            obj = obj.cutScansToSize;
            obj.nScan = cellfun(@numel,obj.TempDataFileObj.TimeCells);
            progressBar.Value = 0.33;

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

            % Average BLK
            if obj.useBlankSubtraction == true && nBlanks > 1
                obj = obj.averageBlankFiles;
                nData = size(obj.TempDataFileObj.ROICells,1); % update number of matrices
            end

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

            if obj.useBlankSubtraction == true % Separate Blank data from Sample data
                tempBLK = obj.TempDataFileObj.ROICells(end,1);
                obj.TempDataFileObj.ROIMatBLK = sparse(tempBLK{:});
                obj.TempDataFileObj.ROICells(end) = [];
                obj.TempDataFileObj.TimeCells(end) = [];
            end

            % subtract blank before IS normalization
            if obj.useBlankSubtraction == true && obj.internalStandardOrder == "BlankIS"
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
            % BLK Subtraction after IS Correction
            if obj.useBlankSubtraction == true && obj.internalStandardOrder == "ISBlank"
                progressBar.Message = "Subtracting Blank";
                roiDataFile = mat2cell(obj.TempDataFileObj.ROIMat,obj.nScanPadded);
                roiDataBlank = obj.TempDataFileObj.ROIMatBLK;
                parfor iFile = 1:size(roiDataFile,1)
                    roiDataFile{iFile,1} = roiDataFile{iFile,1}-padarray(roiDataBlank,size(roiDataFile{iFile,1},1)-size(roiDataBlank,1),0,'post');
                end
                roiDataFile = vertcat(roiDataFile{:});
                roiDataFile = max(roiDataFile,0);
                hasPeak = any(roiDataFile == 0,1);
                obj.TempDataFileObj.ROIMat = roiDataFile(:,hasPeak);
                obj.TempDataFileObj.ROImzVec(:,~hasPeak) = [];
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
            IntegrationData = obj.findPeaks(IDX);
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


            % Build Storage Arrays and filter by number of occurences
            [Output,obj] = obj.buildFeatureArray(IntegrationData);

            %
            % Output = obj.ConfirmSameFeatureByIsotopeDistribution(Output);

            % Occurrence filter
            Output = obj.occurrenceFilterFeatures(Output);

            %gather MS2 spectra
            Output = obj.gatherMS2Spectra(Output);

            %fill remaining fields
            Output = obj.finalizeFeatureOutput(Output);

            % apply scaling
            progressBar.Message = "Apply scaling";
            Output = obj.groupAndSampleScaling(Output);

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
            %delete Temporary file
            delete(obj.TempDataFile)
            obj.TempDataFile = "";
            close(progressBar)
        end

        %% helper functions

        function obj = identifyInternalStandard(obj)
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
            ISData = obj.buildFeatureArray(ISIntegrationData,obj.ISMassFound);

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

        function IntegrationResults = findPeaks(obj,Index)
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
                "spectrumMS2",[], ...
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
            IntegrationResults = obj.finalizeIntegrationOutput(IntegrationResults,timeArray);
        end



        function [Output,obj] = buildFeatureArray(obj,IntegrationResults,varargin)

            %preallocate Output struct
            Output = struct(...
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

            Output.fileNames = obj.fileName;
            Output.groupName = obj.groupName;
            Output.separationType = obj.separationType;
            %check for empty IntegrationResults
            if ~isempty(IntegrationResults)
                %store group infos
                Output.minWidthFiltered = sum(vertcat(IntegrationResults(:).minWidthFiltered));
                Output.maxWidthFiltered = sum(vertcat(IntegrationResults(:).maxWidthFiltered));
                Output.entropyFiltered = sum(vertcat(IntegrationResults(:).entropyFiltered));
                Output.signal2NoiseFiltered = sum(vertcat(IntegrationResults(:).signal2NoiseFiltered));

                %remove unnecessary fields from input struct
                IntegrationResults = rmfield(IntegrationResults,["minWidthFiltered","maxWidthFiltered","entropyFiltered","signal2NoiseFiltered"]);
                nFiles = numel(obj.fileName);

                EmptyStruct = struct(...
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
                timeTolerance = obj.peakTimeTolerance;

                if isscalar(varargin)
                    isISIntegration = true;
                else
                    isISIntegration = false;
                end

                %match features and store in feature struct

                parfor iFeature = 1:length(IntegrationResults)
                    currentFeatureStruct = EmptyStruct;
                    currentFeatureStruct.mass_measured = IntegrationResults(iFeature).mass;
                    currentFeatureStruct.XIC = IntegrationResults(iFeature).XIC;

                    nPeaks = numel(vertcat(IntegrationResults(iFeature).peakLocation{:}));

                    %unpack data
                    peakData = zeros(nPeaks,10);
                    peakData(:,1) = vertcat(IntegrationResults(iFeature).peakLocation{:});
                    peakData(:,2) = vertcat(IntegrationResults(iFeature).peakRetentionTime{:});
                    peakData(:,3) = vertcat(IntegrationResults(iFeature).peakStartLocation{:});
                    peakData(:,4) = vertcat(IntegrationResults(iFeature).peakEndLocation{:});
                    peakData(:,5) = vertcat(IntegrationResults(iFeature).peakHeight{:});
                    peakData(:,6) = vertcat(IntegrationResults(iFeature).peakArea{:});
                    peakData(:,7) = vertcat(IntegrationResults(iFeature).entropy{:});
                    peakData(:,8) = vertcat(IntegrationResults(iFeature).signal2Noise{:});
                    peakData(:,9) = vertcat(IntegrationResults(iFeature).fileID{:});
                    peakData(:,10) = (peakData(:,4)-peakData(:,1))./(peakData(:,1)-peakData(:,3)); %asymmetry factor

                    % preallocate current feature Storage
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

                    %remove empty struct
                    id = isnan([currentFeatureStruct(:).retentionTime])';
                    currentFeatureStruct(id) = [];

                    %store currentFeatureStruct
                    storedFeatures{iFeature,1} = currentFeatureStruct;
                end

                %unzip features
                storedFeatures = vertcat(storedFeatures{:});
                %remove asymmetry field
                storedFeatures = rmfield(storedFeatures,"asymmetry");
                if isISIntegration == false
                    %gather original scans
                    storedFeatures = obj.findOriginalMassScans(storedFeatures);
                end

                Output.feature = storedFeatures;
            else
                Output.feature = IntegrationResults;
            end
        end

        function outputStruct = gatherMS2Spectra(obj,outputStruct)
            %check if MSn data is already loaded
            if isscalar(obj.RawDataFileObj.centroidDataMS2)
                obj = obj.readData(obj.dataFile,obj.separationType);
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

            featureArray = outputStruct.feature;

            parfor iFeature = 1:length(featureArray)
                idMass = [];
                switch MASSUNIT
                    case "Da"
                        idMass = abs(precursor-featureArray(iFeature).mass_measured) <= MASSTOLERANCE;
                    case "ppm"
                        idMass = abs(precursor-featureArray(iFeature).mass_measured)./featureArray(iFeature).mass_measured*10^6 <= MASSTOLERANCE;
                end
                idTime = abs(timeArray-featureArray(iFeature).retentionTime) <= TIMETOLERANCE;
                id = idTime & idMass;
                foundScan = scanArray(id);
                %remove possible empty scans
                foundScan(cellfun(@isempty, foundScan)) = [];
                if numel(foundScan) >= 1
                    foundScan = alignSpectra(foundScan,"average","low","true");
                    foundScan = cleanScans(foundScan,"threshold",0.05);
                else % no found scan
                    foundScan = [];
                end
                featureArray(iFeature).spectrumMS2 = foundScan;
            end
            outputStruct.feature = featureArray;
        end
    end
end