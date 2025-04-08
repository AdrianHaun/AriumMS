classdef GCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        separationType (1,1) string = "GC"
    end

    methods
        function obj = GCData(groupName,window)
            %Construct an instance of this class
            if nargin == 0
                groupName = 0;
                window = 0;
            end
            obj = obj@RawData(groupName,window);
            % set default parameters
            obj.withinFileMassTolerance = 0.1;
            obj.withinFileMassUnit = "Da";
            obj.roiMinOccurence = 10;
            obj.peakMinWidth = 0.8;
            obj.peakMaxWidth = 10;
            obj.betweenFileMassTolerance = 0.05;
            obj.betweenFileMassUnit = "Da";
        end


        %% Data Processing
        function [outputFeatureStruct,obj] = extractFeaturesFromMassData(obj,varargin)

            %check if old results exist and delete them
            if isfile(obj.ROIDataFile)
                delete(obj.ROIDataFile)
            end
            % check for OptimizationMode
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
                fileArray = [fileArray;obj.blankFile];
            end
            %remove possible empty
            iFile = cellfun(@isempty,fileArray);
            fileArray(iFile) = [];
            nData = numel(fileArray);

            progressBar.Message = "Loading files";
            %check if files already loaded then skip loading stage
            test = obj.RawDataFileObj.centroidedDataMS1(1,1);
            if isempty(test{1,1}) || size([obj.fileName;obj.blankFile],1) ~= size(obj.RawDataFileObj.centroidedDataMS1,1)
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
            progressBar.Value = 0.33;
            obj.nScan = cellfun(@numel,obj.TempDataFileObj.TimeCells);
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
            % BLK Subtraction after IS Correction
            if obj.useBlankSubtraction == true && obj.internalStandardOrder == "ISBlank"
                progressBar.Message = "Subtracting Blank";
                roiDataFile = mat2cell(obj.TempDataFileObj.ROIMat,obj.nScansPadded);
                roiDataBlank = obj.TempDataFileObj.ROIMatBLK;
                parfor iFile = 1:size(roiDataFile,1)
                    roiDataFile{iFile,1} = roiDataFile{iFile,1}-padarray(roiDataBlank,size(roiDataFile{iFile,1},1)-size(roiDataBlank,1),0,'post');
                end
                roiDataFile = vertcat(roiDataFile{:});
                roiDataFile = max(roiDataFile,0);
                id = all(roiDataFile >= obj.thresh,1);
                obj.TempDataFileObj.ROIMat = roiDataFile(:,id);
                obj.TempDataFileObj.ROImzVec(:,~id) = [];
                progressBar.Value = progressBar.Value + 0.05;
            end
            % mass correction
            if obj.MassCal == true && ~isempty(obj.ISValue)
                progressBar.Message = "Performing IS mass correction";
                obj = obj.useISMassCorrection;
                progressBar.Value = progressBar.Value + 0.05;
            end

            % Integrate all Peaks
            progressBar.Message = "Picking Peaks";
            IntegrationData = obj.findPeaks;
            progressBar.Value = 0.9;

            progressBar.Message = "Processing found Features";
            IntegrationData = obj.gatherEISpectra(IntegrationData);
            IntegrationData = obj.assignFileID(IntegrationData);
            IntegrationData = obj.fileSortPeaks(IntegrationData);
            IntegrationData = obj.mergeDuplicatePeaksWithinFile(IntegrationData);

            % Build Storage Arrays
            outputFeatureStruct = obj.buildFeatureArray(IntegrationData);

            % confirm same feature by MS2 comparison
            outputFeatureStruct = obj.confirmSameFeatureByMS2(outputFeatureStruct);

            % Occurrence filter
            outputFeatureStruct = obj.occurrenceFilterFeatures(outputFeatureStruct);

            %build average EI (MS2) spectrum
            outputFeatureStruct.feature = obj.finalizeEISpectra(outputFeatureStruct.feature);

            %fill remaining fields
            outputFeatureStruct = obj.finalizeFeatureOutput(outputFeatureStruct);
            progressBar.Value = 0.95;

            %apply scaling
            progressBar.Message = "Apply scaling";
            outputFeatureStruct = obj.groupAndSampleScaling(outputFeatureStruct);

            progressBar.Message = "Group processing successful";
            progressBar.Value = 1;
            
            obj.Output = outputFeatureStruct;

            %post processing cleanup
            if ~exist("mode","var") %save results if batch mode
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

        %% Processing functions
        function IntegrationResults = findPeaks(obj)
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
            %tics = mat2cell(tics,obj.nScansPadded);
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
            IntegrationResults = obj.finalizeIntegrationOutput(IntegrationResults,currentTime);
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

        function output = buildFeatureArray(obj,IntegrationResults,varargin)
            nFiles = numel(obj.fileName);

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

            %store group info
            output.minWidthFiltered = IntegrationResults(1).minWidthFiltered;
            output.maxWidthFiltered = IntegrationResults(1).maxWidthFiltered;
            output.entropyFiltered = IntegrationResults(1).entropyFiltered;
            output.signal2NoiseFiltered = IntegrationResults(1).signal2NoiseFiltered;
            output.fileNames = obj.fileName;
            output.groupName = obj.groupName;
            output.separationType = obj.separationType;

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
                        currentSpectrum = spectrum(~cellfun(@isempty, spectrum));
                        currentSpectrum = vertcat(currentSpectrum{:});
                        %build current average spectrum
                        currentSpectrum = AlignSpectra(currentSpectrum,"average");
                        %gather possible spectra
                        possibleSpectra = IntegrationResults(iFile).spectrumMS2(idx);
                        possibleSpectra = AlignSpectra(possibleSpectra,"normal");
                        %calculate scores
                        [CompoundScores,~] = ScoresBetweenSets(currentSpectrum,possibleSpectra);
                        if all(CompoundScores(:,1) < 700) %features don´t match
                            continue
                        end

                    elseif sum(idx) > 1 %use feat with higher Similarity score
                        currentSpectrum = spectrum(~cellfun(@isempty, spectrum));
                        currentSpectrum = vertcat(currentSpectrum{:});
                        %build current average spectrum
                        currentSpectrum = AlignSpectra(currentSpectrum,"average");
                        %gather possible spectra
                        possibleSpectra = IntegrationResults(iFile).spectrumMS2(idx);
                        possibleSpectra = AlignSpectra(possibleSpectra,"normal");
                        %calculate scores
                        [CompoundScores,~] = ScoresBetweenSets(currentSpectrum,possibleSpectra);
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
                    spectrum{1,iFile} = IntegrationResults(iFile).spectrumMS2(idx);
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

            output.feature = featureStruct;
        end

    end
    methods (Static)

        function outputStruct = finalizeEISpectra(inputStruct)
            %% merges all found EI fragment spectra (all files) into an average spectrum
            outputStruct = inputStruct;
            for n = 1:length(inputStruct)
                spectraCells = inputStruct(n).spectrumMS2;
                spectraCells = horzcat(spectraCells{:});
                outputStruct(n).spectrumMS2 = alignSpectra(spectraCells,"average","low","true");
            end
        end
    end
end