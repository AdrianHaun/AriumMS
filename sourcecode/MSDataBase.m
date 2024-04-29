classdef MSDataBase
    %MSDATABASE Summary of this class goes here
    %   Detailed explanation goes here

    properties

        DataBaseFile    (1,1) string
        Features        (:,:) double
        Significant     (:,1) logical
        HighFoldChange  (:,1) logical
        FoundInGroup    (:,:) string
        IonizationNames (5,1) string
        %UI elements
        Window
        Labels
        SelectDBFile
        GenerateDB
        MS1SearchButton
        MS2SearchButton
        mzTolEditfield
        ErrorUnitDropDown
        FilterFragmentationEnergyCheckbox
        FilterIonizationTypeCheckbox
        FilterSignificantCheckbox
        FilterFoldChangeCheckbox
        DatabaseTypeSwitch
        EICheckBox
        CICheckBox
        ESICheckBox
        APCICheckBox
        APPICheckBox
    end

    methods
        function obj = MSDataBase(CallingApp)
            %MSDATABASE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Features = CallingApp.FeatureData.IdentifierArray;
            if size(CallingApp.FeatureData.GroupName,2) < 2
                obj.Significant = true(size(CallingApp.FeatureData.AverageIntensities));
                obj.HighFoldChange = true(size(CallingApp.FeatureData.AverageIntensities));
            else
                IsSignificant = horzcat(CallingApp.FeatureData.SignificantFeature{:});
                IsSignificant = any(IsSignificant <= CallingApp.maxP,2);
                obj.Significant = IsSignificant;
                Fold = max(horzcat(CallingApp.FeatureData.FullFoldChanges{:}),[],2);
                obj.HighFoldChange = Fold >= CallingApp.minFold;
            end
            obj.FoundInGroup = CallingApp.FeatureData.InGroup;
            obj.IonizationNames = ["EI";"CI";"ESI";"APCI";"APPI"];
            obj.DataBaseFile = CallingApp.dataBaseFile;
            %build UIelements
            obj.Window = uifigure("WindowStyle","normal",...
                "Name","Database Searcher",...
                "Position",[200,200,500,250],...
                "NumberTitle","off");
            obj.Labels(1) = uilabel(obj.Window,"Text","Local Database:","Position",[20,200,160,40],HorizontalAlignment="left",FontName='Impact',FontSize=22);
            obj.SelectDBFile = uibutton(obj.Window,"push",...
                "Position",[170,200,155,40],...
                Text="Select Database File", ...
                Enable="on");
            obj.GenerateDB = uibutton(obj.Window,"push",...
                "Position",[335,200,155,40],...
                Text="Build Database File", ...
                Enable="on");
            obj.MS1SearchButton = uibutton(obj.Window,"push",...
                "Position",[335,90,155,55],...
                Text="MS1 search", ...
                Enable="off");
            obj.MS2SearchButton = uibutton(obj.Window,"push",...
                "Position",[335,25,155,55],...
                Text="MS2 search", ...
                Enable="off");
            obj.mzTolEditfield = uieditfield(obj.Window,"numeric",...
                "Position",[205,145,50,20],...
                "Limits",[0,Inf],...
                "LowerLimitInclusive","off",...
                "Value",0.01,...
                Tooltip="Lowest mass tolerance to match features to database entry.");
            obj.Labels(2) = uilabel(obj.Window,"Text","Allowed mass tolerance:","Position",[60,145,140,20],HorizontalAlignment="left");
            obj.ErrorUnitDropDown = uidropdown(obj.Window, ...
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
            obj.FilterSignificantCheckbox = uicheckbox(obj.Window,...
                "Value",false,...
                "Position",[60,80,250,20],...
                "Text","Search only significant features",...
                Enable="off",...
                Tooltip="Only consider significant features in database search.");
            obj.FilterFoldChangeCheckbox = uicheckbox(obj.Window,...
                "Value",false,...
                "Position",[60,105,250,20],...
                "Text","Search only features with high fold change",...
                Enable="off",...
                Tooltip="Only consider features with fold change greater than minimum fold change in database search.");
            obj.DatabaseTypeSwitch = uiswitch(obj.Window,...
                "Items", {'MassBank','HMDB'},"ItemsData",{0,1},"Position",[405,180,155,15],Enable="on");

            obj.FilterIonizationTypeCheckbox = uicheckbox(obj.Window,...
                "Value",false,...
                "Position",[60,55,250,20],...
                "Text","Ionization Types",...
                Enable="off",...
                Tooltip="Only consider database MS2 spectra with matching ionization type.");

            obj.EICheckBox = uicheckbox(obj.Window,...
                "Value",true,...
                "Position",[90,35,50,20],...
                "Text","EI",...
                Enable="off",...
                Tooltip="Includes EI spectra in Database search results.");
            obj.CICheckBox = uicheckbox(obj.Window,...
                "Value",true,...
                "Position",[90,15,50,20],...
                "Text","CI",...
                Enable="off",...
                Tooltip="Includes CI spectra in Database search results.");
            obj.ESICheckBox = uicheckbox(obj.Window,...
                "Value",true,...
                "Position",[150,35,50,20],...
                "Text","ESI",...
                Enable="off",...
                Tooltip="Includes ESI spectra in Database search results.");
            obj.APCICheckBox = uicheckbox(obj.Window,...
                "Value",true,...
                "Position",[150,15,50,20],...
                "Text","APCI",...
                Enable="off",...
                Tooltip="Includes APCI spectra in Database search results.");
            obj.APPICheckBox = uicheckbox(obj.Window,...
                "Value",true,...
                "Position",[210,35,50,20],...
                "Text","APPI",...
                Enable="off",...
                Tooltip="Includes APPI spectra in Database search results.");

            %add function callbacks
            obj.FilterIonizationTypeCheckbox.ValueChangedFcn = @(src,event) {IonTypeFilterSwitch(obj,src,event)};
            drawnow
            obj = obj.CheckStatus;
        end

        function obj=testDBfile(obj)          %check for existing SQLlite database
            [~,~,ext] = fileparts(obj.DataBaseFile);
            if ~isfile(obj.DataBaseFile) || ~strcmp(ext,".db")
                obj = obj.GenerateDatabase;
            end
        end

        function obj = IonTypeFilterSwitch(obj,~,event)

            if event.Value == true
                %enable sub checkboxes
                obj.EICheckBox.Enable = "on";
                obj.CICheckBox.Enable = "on";
                obj.ESICheckBox.Enable = "on";
                obj.APCICheckBox.Enable = "on";
                obj.APPICheckBox.Enable = "on";
            else
                %disable sub checkboxes
                obj.EICheckBox.Enable = "off";
                obj.CICheckBox.Enable = "off";
                obj.ESICheckBox.Enable = "off";
                obj.APCICheckBox.Enable = "off";
                obj.APPICheckBox.Enable = "off";
            end
        end

        function obj = SelectDB(obj,~,~)
            [file,path] = uigetfile('*.db','Select database file');
            figure(obj.Window);
            %check for user cancel
            if file == 0
                uialert(obj.Window,"Task aborted by user","No database selected")
                obj = obj.CheckStatus;
            else
                obj.DataBaseFile=fullfile(path,file);
                obj = obj.testDBfile;
                obj = obj.CheckStatus;
            end
        end
        function obj = GenerateDatabase(obj,~,~)
            if obj.DataBaseFile == ""
                [file,path] = uiputfile('*.db','Select Save Location and Filename for database file',"DataBase.db");
                figure(obj.Window);
                %check for user cancel
                if file == 0
                    uialert(obj.Window,"Task aborted by user","No database file created","Icon","warning")
                    return
                end
                obj.DataBaseFile=fullfile(path,file);
            else
                CurrentDB = obj.DataBaseFile;
                selection = uiconfirm(obj.Window,"Overwrite "+CurrentDB+" ?","Confirm","Options",["Overwrite","Save as new","Cancel"],"DefaultOption",2,"CancelOption",3);

                switch selection
                    case "Overwrite"
                        %delete file
                        delete(CurrentDB)
                    case "Save as new"
                        [file,path] = uiputfile('*.db','Select Save Location and Filename for database file',"DataBase.db");
                        figure(obj.Window);
                        obj.DataBaseFile=fullfile(path,file);
                    case "Cancel"
                        uialert(obj.Window,"Task aborted by user","No database file created","Icon","warning")
                        return
                end

            end
            DBType = obj.DatabaseTypeSwitch.Value;

            switch DBType
                case 0
                    obj.GenerateMassBankDatabase;
                case 1
                    obj.GenerateHMDBDatabase;
            end

        end

        function DecodedMS2Data = DecodeString(obj,EncodedStrings)
            DecodedMS2Data=cell(size(EncodedStrings,1),1);
            parfor n=1:size(EncodedStrings,1)
                Decoded = matlab.net.base64decode(EncodedStrings(n));
                Decoded = typecast(Decoded,'double');
                DecodedMS2Data{n} = reshape(Decoded,[],2);
            end
        end

        function obj = GenerateMassBankDatabase(obj)
            %clear database if one exists
            if isfile(obj.DataBaseFile)
                delete(obj.DataBaseFile)
            end

            MassBankRoot = uigetdir(pwd,'Open MassBank data main folder');
            figure(obj.Window);
            if isempty(MassBankRoot)
                uialert(obj.Window,"Task aborted by user","No database file created","Icon","warning")
                return
            end
            d =  uiprogressdlg(obj.Window,'Title','Building Database',...
                'Indeterminate','on');
            drawnow

            % Use the 'dir' function to list all files and folders in the root directory
            fileList = dir(fullfile(MassBankRoot, '**', '*.txt'));
            % Create an empty cell array to store your data
            database = struct('ACCESSION', {}, 'NAME', {}, 'FORMULA', {}, ...
                'EXACT_MASS', {}, 'INSTRUMENT_TYPE', {},  'IONIZATION', {}, 'FRAGMENTATION_ENERGY', {},'SPECTRUM', {});
            % Loop through the fileList and read the text files
            parfor i = 1:numel(fileList)
                filename = strcat(fileList(i).folder,'\',fileList(i).name);
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
                for j = 1:numel(fileLines)
                    line = fileLines{j};
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
                        t = j;
                        while t < numel(fileLines)
                            t = t + 1;
                            peakLine = strtrim(fileLines{t});
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
                peaks=reshape(peaks,1,[]);
                peaks=typecast(peaks,'uint8');
                peaks=matlab.net.base64encode(peaks);
                % Add the extracted data to the database structure
                entry = struct('ACCESSION', accession, 'NAME', name, ...
                    'FORMULA', formula, 'EXACT_MASS', exact_mass, ...
                    'INSTRUMENT_TYPE', instrument_type, 'IONIZATION', ion_mode, ...
                    'FRAGMENTATION_ENERGY', FragEnergy, 'SPECTRUM', peaks);

                database = [database, entry];
            end
            database=struct2table(database);
            DataBaseConnection = sqlite(obj.DataBaseFile,"create");
            sqlwrite(DataBaseConnection,"SpectralData",database);
            close(DataBaseConnection)
            close(d)
            obj = obj.CheckStatus;
            uialert(obj.Window,"Local database generation successfull","Database file created","Icon","success")
        end

        function obj = GenerateHMDBDatabase(obj)
            [file,path]=uigetfile('*.xml',"Select HMDB metabolite file (.xml)","MultiSelect","off");

            MetaboliteFile = fullfile(path,file);
            if isempty(MetaboliteFile)
                uialert(obj.Window,"Task aborted by user","No database file created","Icon","warning")
                return
            end

            path = uigetdir(path,"Select HMDB MS/MS Spectral File Folder");
            if isempty(path)
                uialert(obj.Window,"Task aborted by user","No database file created","Icon","warning")
                return
            end

            if exist(obj.DataBaseFile, 'file') == 2
                delete(obj.DataBaseFile);
            end
            d =  uiprogressdlg(obj.Window,'Title','Building Database',...
                'Message',"Loading " + MetaboliteFile,...
                'Indeterminate','on');
            drawnow

            database = struct('ACCESSION', {}, 'NAME', {}, 'FORMULA', {}, ...
                'EXACT_MASS', {}, 'INSTRUMENT_TYPE', {},  ...
                'IONIZATION', {},'FRAGMENTATION_ENERGY', {}, 'SPECTRUM', {});

            %% read Metabolite file
            doc = xmlread(MetaboliteFile);

            MetaboliteList = doc.getElementsByTagName('metabolite');
            numMetabolites = MetaboliteList.getLength;


            %% create file and add first entry
            n=0;
            x = n+1;
            d = uiprogressdlg(obj.Window,'Title','Building Database',...
                'Message',"Metabolite " + x + " of " +numMetabolites);
            d.Value = 0;
            Metabolite = MetaboliteList.item(n);
            %HMDB identifier
            Accession = Metabolite.getElementsByTagName('accession');
            Accession = string(Accession.item(0).getTextContent);

            %Compound Name
            Name = Metabolite.getElementsByTagName('name');
            Name = string(Name.item(0).getTextContent);

            %Compound Mass
            ExactMass = Metabolite.getElementsByTagName('monisotopic_molecular_weight');
            ExactMass = str2double(ExactMass.item(0).getTextContent);

            %Compound Formula
            propertyList = Metabolite.getElementsByTagName('property');
            found = false;
            P = 0;
            Formula = strings(1);
            while found == false & P <= propertyList.getLength-1
                property = propertyList.item(P);
                type = property.getElementsByTagName('kind');
                type = string(type.item(0).getTextContent);
                if strcmp(type,"formula")
                    prop = property.getElementsByTagName('value');
                    Formula = string(prop.item(0).getTextContent);
                    found = true;
                end
                P = P+1;
            end
            % identify correct spectra files
            fstruct = dir(path + "\"+Accession+"*");
            numSpectra = numel(fstruct);
            %load spectra files
            if numSpectra > 0
                parfor S = 1:numSpectra

                    SpectraFile = fullfile(path, fstruct(S).name);
                    [instrumenttype,ionization_mode,collisionenergy,encodedSpectrum] = readHMDBSpectraFile(SpectraFile);
                    entry = struct('ACCESSION', Accession, 'NAME', Name, ...
                        'FORMULA', Formula, 'EXACT_MASS', ExactMass, ...
                        'INSTRUMENT_TYPE', instrumenttype, 'IONIZATION', ionization_mode, ...
                        'FRAGMENTATION_ENERGY', collisionenergy, 'SPECTRUM', encodedSpectrum);
                    database = [database, entry];
                end
            end
            database=struct2table(database);
            DataBaseConnection = sqlite(obj.DataBaseFile,"create");
            sqlwrite(DataBaseConnection,"SpectralData",database);
            database = [];
            d.Value = 1/numMetabolites;
            %% add remaining entries
            for n=1:numMetabolites-1
                x = n+1;
                d.Message = "Metabolite " + x + " of " +numMetabolites;
                Metabolite = MetaboliteList.item(n);
                %HMDB identifier
                Accession = Metabolite.getElementsByTagName('accession');
                Accession = string(Accession.item(0).getTextContent);

                %Compound Name
                Name = Metabolite.getElementsByTagName('name');
                Name = string(Name.item(0).getTextContent);

                %Compound Mass
                ExactMass = Metabolite.getElementsByTagName('monisotopic_molecular_weight');
                ExactMass = str2double(ExactMass.item(0).getTextContent);

                %Compound Formula
                propertyList = Metabolite.getElementsByTagName('property');
                found = false;
                P = 0;
                Formula = strings(1);
                while found == false & P <= propertyList.getLength-1
                    property = propertyList.item(P);
                    type = property.getElementsByTagName('kind');
                    type = string(type.item(0).getTextContent);
                    if strcmp(type,"formula")
                        prop = property.getElementsByTagName('value');
                        Formula = string(prop.item(0).getTextContent);
                        found = true;
                    end
                    P = P+1;
                end
                % identify correct spectra files
                fstruct = dir(path + "\"+Accession+"*");
                numSpectra = numel(fstruct);
                %load spectra files
                if numSpectra > 0
                    parfor S = 1:numSpectra
                        SpectraFile = fullfile(path, fstruct(S).name);
                        [instrumenttype,ionization_mode,collisionenergy,encodedSpectrum] = readHMDBSpectraFile(SpectraFile);
                        entry = struct('ACCESSION', Accession, 'NAME', Name, ...
                            'FORMULA', Formula, 'EXACT_MASS', ExactMass, ...
                            'INSTRUMENT_TYPE', instrumenttype, 'IONIZATION', ionization_mode, ...
                            'FRAGMENTATION_ENERGY', collisionenergy, 'SPECTRUM', string(encodedSpectrum));
                        database = [database, entry];
                    end
                    database=struct2table(database);
                    sqlwrite(DataBaseConnection,"SpectralData",database);
                    database = [];
                    d.Value = (n+1)/numMetabolites;
                    drawnow
                end
            end

            close(DataBaseConnection)
            close(d)
            obj = obj.CheckStatus;
            uialert(obj.Window,"Local database generation successfull","Database file created","Icon","success")

        end
        function out = DataBaseMS1Query(obj)
            QueryMasses = obj.Features(:,1);
            Modifier = true(size(QueryMasses));
            if obj.FilterSignificantCheckbox.Value == true
                Modifier = obj.Significant;
            end
            if obj.FilterFoldChangeCheckbox.Value == true
                Modifier = [Modifier,obj.HighFoldChange];
            end
            Modifier = all(Modifier,2);
            QueryMasses(~Modifier)=[];
            QueryMasses = unique(QueryMasses);
            switch obj.ErrorUnitDropDown.Value
                case "Da"
                    MZmin = QueryMasses - obj.mzTolEditfield.Value;
                    MZmax = QueryMasses + obj.mzTolEditfield.Value;
                case "ppm"
                    MZmin = QueryMasses - (QueryMasses*obj.mzTolEditfield.Value*10^-6);
                    MZmax = QueryMasses + (QueryMasses*obj.mzTolEditfield.Value*10^-6);
            end

            query = ['SELECT NAME, ' ...
                '	FORMULA, ' ...
                '	EXACT_MASS ' ...
                'FROM SpectralData ' ...
                'WHERE EXACT_MASS <= '];
            query = append(query,convertStringsToChars(MZmax + " AND EXACT_MASS >= " + MZmin));
            databasefile = obj.DataBaseFile;
            QueryResults = cell(size(query));
            parfor n=1:length(QueryMasses)
                local_connection = sqlite(databasefile);
                result = fetch(local_connection, query{n});
                %calculate difference in ppm
                close(local_connection);
                if ~isempty(result)
                    result = unique(result,"rows");
                    diff = ((abs(QueryMasses(n)-result.EXACT_MASS))./QueryMasses(n))*10^6;
                    result.DELTA_ppm = round(diff,2);
                    result = sortrows(result,"DELTA_ppm","ascend");
                    QueryResults{n}=result;
                end
            end
            %sort data to features
            out = cell(length(obj.Features(:,1)),1);
            for n=1:size(QueryResults,1)
                idx = obj.Features(:,1) == QueryMasses(n);
                out(idx) = QueryResults(n);
            end
            out(~Modifier) = {};
        end

        function [ResultStorage,DBSpectra,SpectraIndexStorage] = DataBaseMS2Query(obj,MeasuredSpectra)
            %Query time
            databasefile = obj.DataBaseFile;
            QueryMasses = obj.Features(:,1);
            mzTol = obj.mzTolEditfield.Value;
            ErrorUnit = obj.ErrorUnitDropDown.Value;
            Modifier = true(size(QueryMasses));
            if obj.FilterSignificantCheckbox.Value == true
                Modifier = obj.Significant;
            end
            if obj.FilterFoldChangeCheckbox.Value == true
                Modifier = [Modifier,obj.HighFoldChange];
            end
            Modifier = all(Modifier,2);
            NoSpectra = cellfun(@isempty,MeasuredSpectra);
            Modifier(NoSpectra) = false;
            AllowedIonization = obj.IonizationNames([obj.EICheckBox.Value;obj.CICheckBox.Value;obj.ESICheckBox.Value;obj.APCICheckBox.Value;obj.APPICheckBox.Value]);
            switch ErrorUnit
                case "Da"
                    MZmin = QueryMasses - mzTol;
                    MZmax = QueryMasses + mzTol;
                case "ppm"
                    MZmin = QueryMasses - (QueryMasses*mzTol*10^-6);
                    MZmax = QueryMasses + (QueryMasses*mzTol*10^-6);
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
            query = append(query,convertStringsToChars(MZmax + " AND EXACT_MASS >= " + MZmin));
            %databasefile = obj.DataBaseFile;
            QueryResults = cell(size(query));
            local_connection = sqlite(databasefile);
            for n=1:length(QueryMasses)
                if Modifier(n) == true
                    result = fetch(local_connection, query{n});
                else
                    result = [];
                end
                if ~isempty(result)
                    idx = contains(result.INSTRUMENT_TYPE,AllowedIonization);
                    result(~idx,:) = [];
                    diff = ((abs(QueryMasses(n)-result.EXACT_MASS))./QueryMasses(n))*10^6;
                    result.DELTA_ppm = round(diff,2);
                    QueryResults{n}=result;
                end
            end
            close(local_connection);
            %%
            DBSpectra = cell(size(MeasuredSpectra,1),1);

            % prepare Database Spectra
            parfor n=1:length(QueryResults)
                if isempty(QueryResults{n})
                    continue
                else
                    SpectraStrings = QueryResults{n}.SPECTRUM;
                    Spectra = DecodeStrings(SpectraStrings);
                    %rescale intensities
                    Spectra = cellfun(@(x) [x(:,1),x(:,2)./max(x(:,2))],Spectra,'UniformOutput',false);
                    %remove intensities < 5%
                    Spectra = cellfun(@(x) x(x(:,2)>=0.05,:),Spectra,'UniformOutput',false);
                    Spectra = mergeMatricesWithTolerance(Spectra,0.1,"Da");

                    [~,id] = sort(Spectra(:,1),'ascend');
                    DBSpectra{n} = Spectra(id,:);
                end
            end
            %%
            %calculate Scores for each Group
            Scores = cell(1,size(MeasuredSpectra,2));
            parfor n = 1:size(MeasuredSpectra,2)
                Spectra = [MeasuredSpectra(:,n), DBSpectra];
                Scores{1,n} = OuterFeatScores(Spectra,0.015,"Da");
            end
            Scores = horzcat(Scores{:});
            %% get indices of high score spectra
            ResultStorage = cell(size(QueryMasses));
            SpectraIndexStorage = cell(size(QueryMasses));
            for n = 1:length(QueryResults)
                %merge scores
                SC = vertcat(Scores{n,:});
                if isempty(SC) || isempty(QueryResults{n})
                    continue
                else
                    %get Entry index and Score
                    idx = SC(:,1) >= 800;
                    Val = SC(idx,:);
                    %filter duplicates and sort
                    Val = unique(Val,"rows");
                    [~,idx] = sort(Val(:,1),'ascend');
                    Val = Val(idx,:);
                    %sort Query Results
                    Result = QueryResults{n,1}(Val(:,3),:);
                    Result.SPECTRUM = [];
                    Result.ID_MATCH_SCORE = Val(:,1);
                    [Result,id] = unique(Result,"rows");
                    Val = Val(id,:);
                    if ~isempty(Val)
                        [~,id] = sort(Result.ID_MATCH_SCORE,'descend');
                        ResultStorage{n,1} = Result(id,:);
                        SpectraIndexStorage{n,1} = Val(id,:);
                    end
                end
            end

        end

        function obj = CheckStatus(obj)
            if strcmp(obj.DataBaseFile,"")
                obj.MS1SearchButton.Enable = "off";
                obj.MS2SearchButton.Enable = "off";
                obj.mzTolEditfield.Enable = "off";
                obj.ErrorUnitDropDown.Enable = "off";
                % obj.FilterFragmentationTypeCheckbox.Enable = "off";
                %obj.FilterFragmentationEnergyCheckbox.Enable = "off";
                obj.FilterSignificantCheckbox.Enable = "off";
                obj.FilterFoldChangeCheckbox.Enable = "off";
                obj.FilterIonizationTypeCheckbox.Enable = "off";
            else
                obj.MS1SearchButton.Enable = "on";
                obj.MS2SearchButton.Enable = "on";
                obj.mzTolEditfield.Enable = "on";
                obj.ErrorUnitDropDown.Enable = "on";
                % obj.FilterFragmentationTypeCheckbox.Enable = "on";
                %obj.FilterFragmentationEnergyCheckbox.Enable = "on";
                obj.FilterSignificantCheckbox.Enable = "on";
                obj.FilterFoldChangeCheckbox.Enable = "on";
                obj.FilterIonizationTypeCheckbox.Enable = "on";
            end

            drawnow
        end
    end
end

