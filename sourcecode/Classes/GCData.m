classdef GCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        SeparationType (1,1) string = "GC"
    end

    methods
        function obj = GCData(groupNumber,window)
            %Construct an instance of this class
            if nargin == 0
                groupNumber = 0;
                window = 0;
            end
            obj = obj@RawData(groupNumber,window);
            % set default parameters
            obj.mzerror = 0.1;
            obj.mzErrorUnit = "Da";
            obj.minroi = 10;
            obj.minWidth = 0.8;
            obj.maxWidth = 10;
            obj.mzTol = 0.05;
            obj.mzTolUnit = "Da";
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
            test = obj.RawDataFileObj.centroidedDataMS1(1,1);
            if isempty(test{1,1}) || size([obj.Files;obj.BlankFiles],1) ~= size(obj.RawDataFileObj.centroidedDataMS1,1)
                obj = obj.ReadData(FileLocs,obj.SeparationType);
            end
            progressBar.Value = 0.33;
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
            progressBar.Message = "Integrating Peaks";
            IntegrationData = obj.GCIntegrate;
            progressBar.Value = 0.9;

            progressBar.Message = "Processing found Features";
            IntegrationData = obj.GatherEISpectra(IntegrationData);
            IntegrationData = obj.AssignRT2SampleFile(IntegrationData);
            IntegrationData = obj.FileSortPeaks(IntegrationData);
            IntegrationData = obj.mergeDuplicatePeaksWithinFile_GC(IntegrationData);



            % Build Storage Arrays
            Output = obj.BuildStorageArrays_GC(IntegrationData);

            % confirm same feature by MS2 comparison
            Output = obj.ConfirmSameFeatureByMS2(Output);

            % Occurence filter
            Output = obj.OccurenceFilterFeatures(Output);

            %build average EI (MS2) spectrum
            Output.feature = obj.FinalizeEISpectra(Output.feature);

            %fill remaining fields
            Output = obj.FinalizeBatchOutput(Output);
            progressBar.Value = 0.95;

            %apply scaling
            progressBar.Message = "Apply scaling";
            Output = obj.GroupAndSampleScaling(Output);

            progressBar.Message = "Group processing successful";
            progressBar.Value = 1;

            obj.Output = Output;

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
            %delete Temprorary file
            delete(obj.TempDataFile)
            obj.TempDataFile = "";
            close(progressBar)
        end

        %% helper functions
        function IntResults = GCIntegrate(obj)
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
            IntResults.XIC = [times,tics];
            %gather parameters
            currentTIC = full(tics);
            currentTime = full(times);
            smoothedTIC = smoothdata(currentTIC,"gaussian",4,"omitnan");
            %calculate noise level (10th percentile of non zero values)
            noise = prctile(currentTIC(currentTIC > 0),10);
            [~,peakLoc,peakWidth] = findpeaks(currentTIC,"WidthReference","halfheight");
            %calculate initial borders and bring in correct form
            lowerBorders = max(floor(peakLoc-peakWidth/2),1); % limit lower peak border to scan 1
            upperBorders = min(ceil(peakLoc+peakWidth/2),numel(currentTIC)); % limit upper peak border to last scan
            peaks = [peakLoc,lowerBorders,upperBorders];
            peaks = CWTBorderCorrection(peaks,currentTIC,smoothedTIC);
            IntResults.peakLocation = peaks(:,1);
            IntResults.peakStartLocation = peaks(:,2);
            IntResults.peakEndLocation = peaks(:,3);
            IntResults.peakHeight = peaks(:,4);
            IntResults = obj.FilterPeaks(IntResults,noise);
            % entropy calculation
            IntResults = obj.CalculatePeakEntropy(IntResults);
            IntResults = obj.FinalizeIntegrationOutput(IntResults,currentTime);
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
                %normalize Spectras
                Spectras = Spectras./max(Spectras,[],"all");
                %remove comumns with mor than 50% empty
                id = (sum(Spectras~=0)/height(Spectras))<0.5
                Spectras(:,id)=[]

                %mean spectra
                Spectras = mean(Spectras);
                Spectras = [ROImz(~id);full(Spectras)]';
                %remove rows with intensity < 0.01
                Spectras(Spectras(:,2)<0.01,:) = [];
                FoundSpectra{n,1} = Spectras;

                %identify molecular mass
                masses = flip(Spectras(:,1));
                hasMolecularMass = false;
                counter = 0;
                while hasMolecularMass == false & counter < numel(masses)
                    counter = counter + 1;
                    possibleFragment = masses(counter)-EIlosses;
                    hasMolecularMass = any(min(abs(masses-possibleFragment'))<0.1);
                end
                if hasMolecularMass
                    MolecularMass(n,1) = masses(counter,1);
                else
                    MolecularMass(n,1) = Spectras(end,1);
                end
            end
            IntResults.spectrumMS2 = FoundSpectra;
            IntResults.mass = round(MolecularMass,1);
        end


        function outputStruct = FinalizeEISpectra(obj,inputStruct)
            %% merges all found EI fragment spectra (all files) into an average spectrum
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

        function outArray = mergeDuplicatePeaksWithinFile_GC(obj,inArray)

            mzTolerance = obj.mzerror;
            mzerrorUnit = obj.mzErrorUnit;
            rtTol = obj.RTTol;

            s = size(inArray.mass);

            outArray = struct( ...
                "mass",cell(s), ...
                "peakLocation",cell(s), ...
                "peakRetentionTime",cell(s), ...
                "peakStartLocation",cell(s), ...
                "peakEndLocation",cell(s), ...
                "peakHeight",cell(s), ...
                "peakArea",cell(s), ...
                "entropy",cell(s), ...
                "signal2Noise",cell(s), ...
                "minWidthFiltered",inArray.minWidthFiltered, ...
                "maxWidthFiltered",inArray.maxWidthFiltered, ...
                "entropyFiltered",inArray.entropyFiltered, ...
                "signal2NoiseFiltered",inArray.signal2NoiseFiltered, ...
                "spectrumMS2",cell(s), ...
                "XIC",inArray.XIC, ...
                "fileID",cell(s));

            fnames = fieldnames(outArray);
            fnames(any(fnames == ["XIC","minWidthFiltered","maxWidthFiltered","entropyFiltered","signal2NoiseFiltered"],2)) = []; %remove names from list to skip field in assignmelt loop

            tic = full(inArray.XIC(:,1));

            for n = 1:numel(inArray.mass) %sample loop
                featureID = [inArray.mass{1,n},inArray.peakRetentionTime{1,n}];

                counter = 0;
                while ~isempty(featureID)
                    counter = counter + 1;
                    % find all features that match the current signiture
                    %mass tolerance
                    switch mzerrorUnit
                        case "Da"
                            idm = abs(featureID(:,1)-featureID(1,1)) <= mzTolerance;
                        case "ppm"
                            idm = abs(featureID(:,1)-featureID(1,1))./featureID(1,1)*10^6 <= mzTolerance;
                    end
                    %time tolerance
                    idt = abs(featureID(:,2)-featureID(1,2)) <= rtTol;
                    id = idm & idt;

                    %% multiple peaks found
                    % remove peaks with peakheight < 3x baseline
                    if sum(id) > 1
                        bordersStart = inArray.peakStartLocation{1,n};
                        bordersEnd = inArray.peakEndLocation{1,n};
                        %remove less prominent peak
                        idx = zeros(size(id));
                        for h = 1:numel(id)
                            if id(h) == false
                                continue
                            else
                                temptic = tic(bordersStart(h):bordersEnd(h),:);
                                idx(h) =  mean([temptic(1);temptic(end)]) / max(temptic);
                            end
                        end
                        %remove feat and return while loop
                        id = id & idx(idx==max(idx));
                        for f = 1:numel(fnames) % loop over each field name
                            inArray.(fnames{f}){1,n}(id) =  [];
                        end

                        %% only one peak remaining
                        % store in output and remove from input
                    elseif sum(id) == 1
                        for f = 1:numel(fnames) % loop over each field name
                            outArray(n).(fnames{f})= vertcat(outArray(n).(fnames{f}),inArray.(fnames{f}){1,n}(id));
                            inArray.(fnames{f}){1,n}(id) =  [];
                        end
                    else %no matching peak

                    end
                    % remove entries from current list
                    featureID(id,:) = [];
                end
            end
        end

        function output = BuildStorageArrays_GC(obj,IntegrationResults,varargin)
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
            output.groupName = obj.GroupName;
            output.separationType = obj.SeparationType;

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
                isISIntegration = true;
            else
                isISIntegration = false;
            end

            %match features and store in feature struct

            uniqueFeatures = [vertcat(IntegrationResults(:).mass),vertcat(IntegrationResults(:).peakRetentionTime)];

            % over preallocat feature Storage
            featureStruct = repmat(featureStruct,height(uniqueFeatures),1);
            n = 0;
            while ~isempty(uniqueFeatures)
                n = n + 1;
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

                    elseif sum(idx) == 1 && ~all(isnan(peakLocation)) %compare spectra if already peak stored in feature
                        currentSpectrum = spectrum(~cellfun(@isempty, spectrum));
                        currentSpectrum = vertcat(currentSpectrum{:});
                        %build current average spectrum
                        currentSpectrum = AlignSpectra(currentSpectrum,"average");
                        %gather possible spectra
                        possibleSpectra = IntegrationResults(file).spectrumMS2(idx);
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
                        possibleSpectra = IntegrationResults(file).spectrumMS2(idx);
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
                    areas(1,file) = IntegrationResults(file).peakArea(idx);
                    heights(1,file) = IntegrationResults(file).peakHeight(idx);
                    retentionTimes(1,file) = IntegrationResults(file).peakRetentionTime(idx);
                    peakLocation(1,file) = IntegrationResults(file).peakLocation(idx);
                    peakBorders(1,file) = IntegrationResults(file).peakStartLocation(idx);
                    peakBorders(2,file) = IntegrationResults(file).peakEndLocation(idx);
                    signal2Noise(1,file) = IntegrationResults(file).signal2Noise(idx);
                    entropy(1,file) = IntegrationResults(file).entropy(idx);
                    spectrum{1,file} = IntegrationResults(file).spectrumMS2(idx);
                    xic{1,file} = full(XIC(peakBorders(1,file):peakBorders(2,file),:));

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

            if isISIntegration == false

                %gather original scans
                featureStruct = obj.FindOriginalScans(featureStruct);
            end

            output.feature = featureStruct;
        end

    end
end