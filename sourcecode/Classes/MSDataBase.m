classdef MSDataBase
    %MSDATABASE Class for handling local MS Databese search

    properties
        dataBaseFile            (1,1) string
        features                (:,:) double
        significantFeatureArray (:,1) logical
        highFoldChangeArray     (:,1) logical
        foundInGroup            (:,:) string
        ionizationNames         (5,1) string
        %UI elements
        window
        labels
        selectDataBaseFile
        generateDataBase
        ms1SearchButton
        ms2SearchButton
        mzTolEditfield
        errorUnitDropDown
        filterFragmentationEnergyCheckbox
        filterIonizationTypeCheckbox
        filterSignificantCheckbox
        filterFoldChangeCheckbox
        databaseTypeSwitch
        electronIonizationCheckBox
        chemicalIonizationCheckBox
        electroSprayIonizationCheckBox
        atmosphericPressureChemicalIonizationCheckBox
        atmosphericPressurePhotoIonizationCheckBox
    end

    methods
        function obj = MSDataBase(CallingApp)

            if nargin > 0

                obj.features = CallingApp.FeatureData.IdentifierArray;
                if size(CallingApp.FeatureData.GroupName,2) < 2
                    obj.significantFeatureArray = true(size(CallingApp.FeatureData.AverageIntensities));
                    obj.highFoldChangeArray = true(size(CallingApp.FeatureData.AverageIntensities));
                else
                    IsSignificant = horzcat(CallingApp.FeatureData.SignificantFeature{:});
                    IsSignificant = any(IsSignificant <= CallingApp.maxP,2);
                    obj.significantFeatureArray = IsSignificant;
                    Fold = max(horzcat(CallingApp.FeatureData.FullFoldChanges{:}),[],2);
                    obj.highFoldChangeArray = Fold >= CallingApp.minFold;
                end
                obj.foundInGroup = CallingApp.FeatureData.InGroup;
                obj.ionizationNames = ["EI";"CI";"ESI";"APCI";"APPI"];
                obj.dataBaseFile = CallingApp.dataBaseFile;
                %build UIelements
                obj.window = uifigure("WindowStyle","normal",...
                    "Name","Database Searcher",...
                    "Position",[200,200,500,250],...
                    "NumberTitle","off");
                obj.labels(1) = uilabel(obj.window,"Text","Local Database:","Position",[20,200,160,40],HorizontalAlignment="left",FontName='Impact',FontSize=22);
                obj.selectDataBaseFile = uibutton(obj.window,"push",...
                    "Position",[170,200,155,40],...
                    Text="Select Database File", ...
                    Enable="on");
                obj.generateDataBase = uibutton(obj.window,"push",...
                    "Position",[335,200,155,40],...
                    Text="Build Database File", ...
                    Enable="on");
                obj.ms1SearchButton = uibutton(obj.window,"push",...
                    "Position",[335,90,155,55],...
                    Text="MS1 search", ...
                    Enable="off");
                obj.ms2SearchButton = uibutton(obj.window,"push",...
                    "Position",[335,25,155,55],...
                    Text="MS2 search", ...
                    Enable="off");
                obj.mzTolEditfield = uieditfield(obj.window,"numeric",...
                    "Position",[205,145,50,20],...
                    "Limits",[0,Inf],...
                    "LowerLimitInclusive","off",...
                    "Value",0.01,...
                    Tooltip="Lowest mass tolerance to match features to database entry.");
                obj.labels(2) = uilabel(obj.window,"Text","Allowed mass tolerance:","Position",[60,145,140,20],HorizontalAlignment="left");
                obj.errorUnitDropDown = uidropdown(obj.window, ...
                    "Items",["ppm","Da"],...
                    ItemsData=["ppm","Da"], ...
                    Value="Da",...
                    Position=[260,145,60,20],...
                    Tooltip="Specify whether a fixed mz error in Da or a relative error in ppm is used.");
                % obj.FilterFragmentationEnergyCheckbox = uicheckbox(obj.Window,...
                %     "Value",false,...
                %     "Position",[60,55,250,20],...
                %     "Text","Filter fragmentation energy",...
                %     Tooltip="Only consider database MS2 spectra with matching fragmentation energy.");
                obj.filterSignificantCheckbox = uicheckbox(obj.window,...
                    "Value",false,...
                    "Position",[60,80,250,20],...
                    "Text","Search only significant features",...
                    Enable="off",...
                    Tooltip="Only consider significant features in database search.");
                obj.filterFoldChangeCheckbox = uicheckbox(obj.window,...
                    "Value",false,...
                    "Position",[60,105,250,20],...
                    "Text","Search only features with high fold change",...
                    Enable="off",...
                    Tooltip="Only consider features with fold change greater than minimum fold change in database search.");
                obj.databaseTypeSwitch = uiswitch(obj.window,...
                    "Items", {'MassBank','HMDB'},"ItemsData",{0,1},"Position",[405,180,155,15],Enable="on");

                obj.filterIonizationTypeCheckbox = uicheckbox(obj.window,...
                    "Value",false,...
                    "Position",[60,55,250,20],...
                    "Text","Ionization Types",...
                    Enable="off",...
                    Tooltip="Only consider database MS2 spectra with matching ionization type.");

                obj.electronIonizationCheckBox = uicheckbox(obj.window,...
                    "Value",true,...
                    "Position",[90,35,50,20],...
                    "Text","EI",...
                    Enable="off",...
                    Tooltip="Includes EI spectra in Database search results.");
                obj.chemicalIonizationCheckBox = uicheckbox(obj.window,...
                    "Value",true,...
                    "Position",[90,15,50,20],...
                    "Text","CI",...
                    Enable="off",...
                    Tooltip="Includes CI spectra in Database search results.");
                obj.electroSprayIonizationCheckBox = uicheckbox(obj.window,...
                    "Value",true,...
                    "Position",[150,35,50,20],...
                    "Text","ESI",...
                    Enable="off",...
                    Tooltip="Includes ESI spectra in Database search results.");
                obj.atmosphericPressureChemicalIonizationCheckBox = uicheckbox(obj.window,...
                    "Value",true,...
                    "Position",[150,15,50,20],...
                    "Text","APCI",...
                    Enable="off",...
                    Tooltip="Includes APCI spectra in Database search results.");
                obj.atmosphericPressurePhotoIonizationCheckBox = uicheckbox(obj.window,...
                    "Value",true,...
                    "Position",[210,35,50,20],...
                    "Text","APPI",...
                    Enable="off",...
                    Tooltip="Includes APPI spectra in Database search results.");

                %add function callbacks
                obj.filterIonizationTypeCheckbox.ValueChangedFcn = @(src,event) {ionTypeFilterSwitch(obj,src,event)};
                drawnow
                obj = obj.checkStatus;
            end
        end

        function obj = testDBfile(obj)          %check for existing SQLlite database
            [~,~,ext] = fileparts(obj.dataBaseFile);
            if ~isfile(obj.dataBaseFile) || ~strcmp(ext,".db")
                obj = obj.generateDatabase;
            end
        end

        function obj = ionTypeFilterSwitch(obj,~,event)

            if event.Value == true
                %enable sub checkboxes
                obj.electronIonizationCheckBox.Enable = "on";
                obj.chemicalIonizationCheckBox.Enable = "on";
                obj.electroSprayIonizationCheckBox.Enable = "on";
                obj.atmosphericPressureChemicalIonizationCheckBox.Enable = "on";
                obj.atmosphericPressurePhotoIonizationCheckBox.Enable = "on";
            else
                %disable sub checkboxes
                obj.electronIonizationCheckBox.Enable = "off";
                obj.chemicalIonizationCheckBox.Enable = "off";
                obj.electroSprayIonizationCheckBox.Enable = "off";
                obj.atmosphericPressureChemicalIonizationCheckBox.Enable = "off";
                obj.atmosphericPressurePhotoIonizationCheckBox.Enable = "off";
            end
        end

        function obj = selectDB(obj,~,~)
            [file,path] = uigetfile('*.db','Select database file');
            figure(obj.window);
            %check for user cancel
            if file == 0
                uialert(obj.window,"Task aborted by user","No database selected")
                obj = obj.checkStatus;
            else
                obj.dataBaseFile=fullfile(path,file);
                obj = obj.testDBfile;
                obj = obj.checkStatus;
            end
        end
        function obj = generateDatabase(obj,~,~)
            if obj.dataBaseFile == ""
                [file,path] = uiputfile('*.db','Select Save Location and Filename for database file',"DataBase.db");
                figure(obj.window);
                %check for user cancel
                if file == 0
                    uialert(obj.window,"Task aborted by user","No database file created","Icon","warning")
                    return
                end
                obj.dataBaseFile = fullfile(path,file);
            else
                currentDB = obj.dataBaseFile;
                selection = uiconfirm(obj.window,"Overwrite "+currentDB+" ?","Confirm","Options",["Overwrite","Save as new","Cancel"],"DefaultOption",2,"CancelOption",3);

                switch selection
                    case "Overwrite"
                        %delete file
                        delete(currentDB)
                    case "Save as new"
                        [file,path] = uiputfile('*.db','Select Save Location and Filename for database file',"DataBase.db");
                        figure(obj.window);
                        obj.dataBaseFile=fullfile(path,file);
                    case "Cancel"
                        uialert(obj.window,"Task aborted by user","No database file created","Icon","warning")
                        return
                end

            end
            DBType = obj.databaseTypeSwitch.Value;

            switch DBType
                case 0
                    obj.generateMassBankDatabase;
                case 1
                    obj.generateHMDBDatabase;
            end

        end

        function decodedMS2Data = decodeString(obj,encodedStrings)
            decodedMS2Data = cell(size(encodedStrings,1),1);
            parfor iString = 1:size(encodedStrings,1)
                Decoded = matlab.net.base64decode(encodedStrings(iString));
                Decoded = typecast(Decoded,'double');
                decodedMS2Data{iString} = reshape(Decoded,[],2);
            end
        end

        function obj = generateMassBankDatabase(obj)
            %clear database if one exists
            if isfile(obj.dataBaseFile)
                delete(obj.dataBaseFile)
            end

            massBankRoot = uigetdir(pwd,'Open MassBank data main folder');
            figure(obj.window);
            if isempty(massBankRoot)
                uialert(obj.window,"Task aborted by user","No database file created","Icon","warning")
                return
            end
            d =  uiprogressdlg(obj.window,'Title','Building Database',...
                'Indeterminate','on');
            drawnow

            % Use the 'dir' function to list all files and folders in the root directory
            fileList = dir(fullfile(massBankRoot, '**', '*.txt'));

            % Create an empty cell array to store data
            database = struct('ACCESSION', {}, 'NAME', {}, 'FORMULA', {}, ...
                'EXACT_MASS', {}, 'INSTRUMENT_TYPE', {},  'IONIZATION', {}, 'FRAGMENTATION_ENERGY', {},'SPECTRUM', {});
            % Loop through the fileList and read the text files
            parfor iFile = 1:numel(fileList)
                filename = strcat(fileList(iFile).folder,'\',fileList(iFile).name);
                % Read the file into a cell array
                fileLines = importdata(filename, '\n');
                % Initialize variables to store information
                accession = '';
                name = '';
                formula = '';
                exact_mass = [];
                instrument_type = '';
                ion_mode = '';
                FragEnergy = [];
                peaks = [];

                % Loop through each line in the file
                for jLine = 1:numel(fileLines)
                    line = fileLines{jLine};
                    chNameCaptured = false;
                    % Extract the information you need from each line
                    if startsWith(line, 'ACCESSION:')
                        accession = strtrim(strrep(line, 'ACCESSION:', ''));
                    elseif ~chNameCaptured && startsWith(line, 'CH$NAME:')
                        name = strtrim(strrep(line, 'CH$NAME:', ''));
                        chNameCaptured = true;  % Set the flag to true to indicate that first CH$NAME has been captured
                    elseif startsWith(line, 'CH$FORMULA:')
                        formula = strtrim(strrep(line, 'CH$FORMULA:', ''));
                    elseif startsWith(line, 'CH$EXACT_MASS:')
                        exact_mass = str2double(strtrim(strrep(line, 'CH$EXACT_MASS:', '')));
                    elseif startsWith(line, 'AC$INSTRUMENT_TYPE:')
                        instrument_type = strtrim(strrep(line, 'AC$INSTRUMENT_TYPE:', ''));
                    elseif startsWith(line, 'AC$MASS_SPECTROMETRY: ION_MODE')
                        ion_mode = strtrim(strrep(line, 'AC$MASS_SPECTROMETRY: ION_MODE', ''));
                    elseif startsWith(line, 'AC$MASS_SPECTROMETRY: COLLISION_ENERGY')
                        FragEnergy = str2double(strtrim(strrep(line, 'AC$MASS_SPECTROMETRY: COLLISION_ENERGY', '')));
                    elseif startsWith(line, 'PK$PEAK:')
                        % Extract and parse peak values from subsequent lines
                        peaks = [];
                        kLine = jLine;
                        while kLine < numel(fileLines)
                            kLine = kLine + 1;
                            peakLine = strtrim(fileLines{kLine});
                            if isempty(peakLine)
                                break;
                            end
                            parts = split(peakLine, ' ');
                            if numel(parts) == 3
                                peaks = [peaks; str2double(parts(1)), str2double(parts(3))];
                            end
                        end
                    end
                end
                %check inputs
                if isempty(ion_mode)
                    ion_mode = NaN;
                end
                if isempty(exact_mass)
                    exact_mass = NaN;
                end
                if isempty(FragEnergy) | isnan(FragEnergy)
                    FragEnergy = 0;
                end
                if isempty(peaks)
                    peaks = [0 0];
                end
                %decode peaks to string
                peaks = reshape(peaks,1,[]);
                peaks = typecast(peaks,'uint8');
                peaks = matlab.net.base64encode(peaks);
                % Add the extracted data to the database structure
                entry = struct('ACCESSION', accession, 'NAME', name, ...
                    'FORMULA', formula, 'EXACT_MASS', exact_mass, ...
                    'INSTRUMENT_TYPE', instrument_type, 'IONIZATION', ion_mode, ...
                    'FRAGMENTATION_ENERGY', FragEnergy, 'SPECTRUM', peaks);

                database = [database, entry];
            end
            database = struct2table(database);
            dataBaseConnection = sqlite(obj.dataBaseFile,"create");
            sqlwrite(dataBaseConnection,"SpectralData",database);
            close(dataBaseConnection)
            close(d)
            obj = obj.checkStatus;
            uialert(obj.window,"Local database generation successfull","Database file created","Icon","success")
        end

        function obj = generateHMDBDatabase(obj)
            [file,path] = uigetfile('*.xml',"Select HMDB metabolite file (.xml)","MultiSelect","off");

            metaboliteFile = fullfile(path,file);
            if isempty(metaboliteFile)
                uialert(obj.window,"Task aborted by user","No database file created","Icon","warning")
                return
            end

            path = uigetdir(path,"Select HMDB MS/MS Spectral File Folder");
            if isempty(path)
                uialert(obj.window,"Task aborted by user","No database file created","Icon","warning")
                return
            end

            if exist(obj.dataBaseFile, 'file') == 2
                delete(obj.dataBaseFile);
            end
            d =  uiprogressdlg(obj.window,'Title','Building Database',...
                'Message',"Loading " + metaboliteFile,...
                'Indeterminate','on');
            drawnow

            database = struct('ACCESSION', {}, 'NAME', {}, 'FORMULA', {}, ...
                'EXACT_MASS', {}, 'INSTRUMENT_TYPE', {},  ...
                'IONIZATION', {},'FRAGMENTATION_ENERGY', {}, 'SPECTRUM', {});

            %% read Metabolite file
            doc = xmlread(metaboliteFile);

            metaboliteList = doc.getElementsByTagName('metabolite');
            nMetabolites = metaboliteList.getLength;


            %% create file and add first entry
            iMetabolite = 0;
            currentMetaboliteNumber = iMetabolite+1;
            d = uiprogressdlg(obj.window,'Title','Building Database',...
                'Message',"Metabolite " + currentMetaboliteNumber + " of " +nMetabolites);
            d.Value = 0;
            metabolite = metaboliteList.item(iMetabolite);
            %HMDB identifier
            accession = metabolite.getElementsByTagName('accession');
            accession = string(accession.item(0).getTextContent);

            %Compound Name
            name = metabolite.getElementsByTagName('name');
            name = string(name.item(0).getTextContent);

            %Compound Mass
            exactMass = metabolite.getElementsByTagName('monisotopic_molecular_weight');
            exactMass = str2double(exactMass.item(0).getTextContent);

            %Compound Formula
            propertyList = metabolite.getElementsByTagName('property');
            found = false;
            iProperty = 0;
            formula = strings(1);
            while found == false & iProperty <= propertyList.getLength-1
                property = propertyList.item(iProperty);
                type = property.getElementsByTagName('kind');
                type = string(type.item(0).getTextContent);
                if strcmp(type,"formula")
                    prop = property.getElementsByTagName('value');
                    formula = string(prop.item(0).getTextContent);
                    found = true;
                end
                iProperty = iProperty+1;
            end
            % identify correct spectra files
            fstruct = dir(path + "\"+accession+"*");
            nSpectra = numel(fstruct);
            %load spectra files
            if nSpectra > 0
                parfor iSpectrum = 1:nSpectra

                    spectraFile = fullfile(path, fstruct(iSpectrum).name);
                    [instrumentType,ionizationMode,collisionEnergy,encodedSpectrum] = readHMDBSpectraFile(spectraFile);
                    entry = struct('ACCESSION', accession, 'NAME', name, ...
                        'FORMULA', formula, 'EXACT_MASS', exactMass, ...
                        'INSTRUMENT_TYPE', instrumentType, 'IONIZATION', ionizationMode, ...
                        'FRAGMENTATION_ENERGY', collisionEnergy, 'SPECTRUM', encodedSpectrum);
                    database = [database, entry];
                end
            end
            database = struct2table(database);
            DataBaseConnection = sqlite(obj.dataBaseFile,"create");
            sqlwrite(DataBaseConnection,"SpectralData",database);
            database = [];
            d.Value = 1/nMetabolites;
            %% add remaining entries
            for iMetabolite = 1:nMetabolites-1
                currentMetaboliteIndex = iMetabolite+1;
                d.Message = "Metabolite " + currentMetaboliteIndex + " of " +nMetabolites;
                metabolite = metaboliteList.item(iMetabolite);
                %HMDB identifier
                accession = metabolite.getElementsByTagName('accession');
                accession = string(accession.item(0).getTextContent);

                %Compound Name
                name = metabolite.getElementsByTagName('name');
                name = string(name.item(0).getTextContent);

                %Compound Mass
                exactMass = metabolite.getElementsByTagName('monisotopic_molecular_weight');
                exactMass = str2double(exactMass.item(0).getTextContent);

                %Compound Formula
                propertyList = metabolite.getElementsByTagName('property');
                found = false;
                iProperty = 0;
                formula = strings(1);
                while found == false & iProperty <= propertyList.getLength-1
                    property = propertyList.item(iProperty);
                    type = property.getElementsByTagName('kind');
                    type = string(type.item(0).getTextContent);
                    if strcmp(type,"formula")
                        prop = property.getElementsByTagName('value');
                        formula = string(prop.item(0).getTextContent);
                        found = true;
                    end
                    iProperty = iProperty+1;
                end
                % identify correct spectra files
                fstruct = dir(path + "\"+accession+"*");
                numSpectra = numel(fstruct);
                %load spectra files
                if numSpectra > 0
                    parfor iSpectrum = 1:numSpectra
                        spectraFile = fullfile(path, fstruct(iSpectrum).name);
                        [instrumentType,ionizationMode,collisionEnergy,encodedSpectrum] = readHMDBSpectraFile(spectraFile);
                        entry = struct('ACCESSION', accession, 'NAME', name, ...
                            'FORMULA', formula, 'EXACT_MASS', exactMass, ...
                            'INSTRUMENT_TYPE', instrumentType, 'IONIZATION', ionizationMode, ...
                            'FRAGMENTATION_ENERGY', collisionEnergy, 'SPECTRUM', string(encodedSpectrum));
                        database = [database, entry];
                    end
                    database = struct2table(database);
                    sqlwrite(DataBaseConnection,"SpectralData",database);
                    database = [];
                    d.Value = (iMetabolite+1)/nMetabolites;
                    drawnow
                end
            end

            close(DataBaseConnection)
            close(d)
            obj = obj.checkStatus;
            uialert(obj.window,"Local database generation successfull","Database file created","Icon","success")

        end

        function out = dataBaseMS1Query(obj)
            queryMasses = obj.features(:,1);
            modifier = true(size(queryMasses));
            if obj.filterSignificantCheckbox.Value == true
                modifier = obj.significantFeatureArray;
            end
            if obj.filterFoldChangeCheckbox.Value == true
                modifier = [modifier,obj.highFoldChangeArray];
            end
            modifier = all(modifier,2);
            queryMasses(~modifier)=[];
            queryMasses = unique(queryMasses);
            switch obj.errorUnitDropDown.Value
                case "Da"
                    mzMin = queryMasses - obj.mzTolEditfield.Value;
                    mzMax = queryMasses + obj.mzTolEditfield.Value;
                case "ppm"
                    mzMin = queryMasses - (queryMasses*obj.mzTolEditfield.Value*10^-6);
                    mzMax = queryMasses + (queryMasses*obj.mzTolEditfield.Value*10^-6);
            end

            query = ['SELECT NAME, ' ...
                '	FORMULA, ' ...
                '	EXACT_MASS ' ...
                'FROM SpectralData ' ...
                'WHERE EXACT_MASS <= '];
            query = append(query,convertStringsToChars(mzMax + " AND EXACT_MASS >= " + mzMin));
            databasefile = obj.dataBaseFile;
            queryResults = cell(size(query));
            parfor iQuery = 1:length(queryMasses)
                local_connection = sqlite(databasefile);
                result = fetch(local_connection, query{iQuery});
                %calculate difference in ppm
                close(local_connection);
                if ~isempty(result)
                    result = unique(result,"rows");
                    diff = ((abs(queryMasses(iQuery)-result.EXACT_MASS))./queryMasses(iQuery))*10^6;
                    result.DELTA_ppm = round(diff,2);
                    result = sortrows(result,"DELTA_ppm","ascend");
                    queryResults{iQuery} = result;
                end
            end
            %sort data to features
            out = cell(length(obj.features(:,1)),1);
            for iQuery=1:size(queryResults,1)
                idx = obj.features(:,1) == queryMasses(iQuery);
                out(idx) = queryResults(iQuery);
            end
            out(~modifier) = {};
        end

        function [resultStorage,dataBaseSpectra,spectraIndexStorage] = dataBaseMS2Query(obj,measuredSpectra)
            %Query time
            databasefile = obj.dataBaseFile;
            queryMasses = obj.features(:,1);
            mzTolerance = obj.mzTolEditfield.Value;
            errorUnit = obj.errorUnitDropDown.Value;
            modifier = true(size(queryMasses));
            if obj.filterSignificantCheckbox.Value == true
                modifier = obj.significantFeatureArray;
            end
            if obj.filterFoldChangeCheckbox.Value == true
                modifier = [modifier,obj.highFoldChangeArray];
            end
            modifier = all(modifier,2);
            noSpectra = cellfun(@isempty,measuredSpectra);
            modifier(noSpectra) = false;
            allowedIonization = obj.ionizationNames([obj.electronIonizationCheckBox.Value;obj.chemicalIonizationCheckBox.Value;obj.electroSprayIonizationCheckBox.Value;obj.atmosphericPressureChemicalIonizationCheckBox.Value;obj.atmosphericPressurePhotoIonizationCheckBox.Value]);
            switch errorUnit
                case "Da"
                    mzMin = queryMasses - mzTolerance;
                    mzMax = queryMasses + mzTolerance;
                case "ppm"
                    mzMin = queryMasses - (queryMasses*mzTolerance*10^-6);
                    mzMax = queryMasses + (queryMasses*mzTolerance*10^-6);
            end
            query = ['SELECT NAME, ' ...
                '	FORMULA, ' ...
                '	EXACT_MASS, ' ...
                '	INSTRUMENT_TYPE, ' ...
                '	IONIZATION, ' ...
                '	FRAGMENTATION_ENERGY, ' ...
                '	SPECTRUM ' ...
                'FROM SpectralData '...
                'WHERE EXACT_MASS <= '];
            query = append(query,convertStringsToChars(mzMax + " AND EXACT_MASS >= " + mzMin));
            %databasefile = obj.DataBaseFile;
            queryResults = cell(size(query));

            parfor iQuery = 1:length(queryMasses)
                local_connection = sqlite(databasefile);
                if modifier(iQuery) == true
                    result = fetch(local_connection, query{iQuery});
                else
                    result = [];
                end
                close(local_connection);
                if ~isempty(result)
                    idx = contains(result.INSTRUMENT_TYPE,allowedIonization);
                    result(~idx,:) = [];
                    diff = ((abs(queryMasses(iQuery)-result.EXACT_MASS))./queryMasses(iQuery))*10^6;
                    result.DELTA_ppm = round(diff,2);
                    queryResults{iQuery} = result;
                end
            end

            %%
            dataBaseSpectra = cell(size(measuredSpectra,1),1);

            % prepare Database Spectra
            parfor iQuery = 1:length(queryResults)
                if isempty(queryResults{iQuery})
                    continue
                else
                    spectraStrings = queryResults{iQuery}.SPECTRUM;
                    spectra = decodeSpectra(spectraStrings);
                    %rescale intensities
                    spectra = cellfun(@(x) [x(:,1),x(:,2)./max(x(:,2))],spectra,'UniformOutput',false);
                    %remove intensities < 5%
                    spectra = cellfun(@(x) x(x(:,2)>=0.05,:),spectra,'UniformOutput',false);
                    spectra = mergeMatricesWithTolerance(spectra,0.1,"Da");

                    [~,id] = sort(spectra(:,1),'ascend');
                    dataBaseSpectra{iQuery} = spectra(id,:);
                end
            end
            %%
            %calculate Scores for each Group
            scores = cell(1,size(measuredSpectra,2));
            parfor iQuery = 1:size(measuredSpectra,2)
                spectra = [measuredSpectra(:,iQuery), dataBaseSpectra];
                scores{1,iQuery} = OuterFeatScores(spectra,0.015,"Da");
            end
            scores = horzcat(scores{:});
            %% get indices of high score spectra
            resultStorage = cell(size(queryMasses));
            spectraIndexStorage = cell(size(queryMasses));
            parfor iQuery = 1:length(queryResults)
                %merge scores
                SC = vertcat(scores{iQuery,:});
                if isempty(SC) || isempty(queryResults{iQuery})
                    continue
                else
                    %get Entry index and Score
                    idx = SC(:,1) >= 800;
                    val = SC(idx,:);
                    if isempty(val)
                        continue
                    end
                    %filter duplicates and sort
                    val = unique(val,"rows");
                    [~,idx] = sort(val(:,1),'ascend');
                    val = val(idx,:);
                    %sort Query Results
                    result = queryResults{iQuery,1}(val(:,2),:);
                    result.SPECTRUM = [];
                    result.ID_MATCH_SCORE = val(:,1);
                    [result,id] = unique(result,"rows");
                    val = val(id,:);
                    if ~isempty(val)
                        [~,id] = sort(result.ID_MATCH_SCORE,'descend');
                        resultStorage{iQuery,1} = result(id,:);
                        spectraIndexStorage{iQuery,1} = val(id,:);
                    end
                end
            end

        end

        function obj = checkStatus(obj)
            if strcmp(obj.dataBaseFile,"")
                obj.ms1SearchButton.Enable = "off";
                obj.ms2SearchButton.Enable = "off";
                obj.mzTolEditfield.Enable = "off";
                obj.errorUnitDropDown.Enable = "off";
                % obj.FilterFragmentationTypeCheckbox.Enable = "off";
                %obj.FilterFragmentationEnergyCheckbox.Enable = "off";
                obj.filterSignificantCheckbox.Enable = "off";
                obj.filterFoldChangeCheckbox.Enable = "off";
                obj.filterIonizationTypeCheckbox.Enable = "off";
            else
                obj.ms1SearchButton.Enable = "on";
                obj.ms2SearchButton.Enable = "on";
                obj.mzTolEditfield.Enable = "on";
                obj.errorUnitDropDown.Enable = "on";
                % obj.FilterFragmentationTypeCheckbox.Enable = "on";
                %obj.FilterFragmentationEnergyCheckbox.Enable = "on";
                obj.filterSignificantCheckbox.Enable = "on";
                obj.filterFoldChangeCheckbox.Enable = "on";
                obj.filterIonizationTypeCheckbox.Enable = "on";
            end

            drawnow
        end
    end
end

