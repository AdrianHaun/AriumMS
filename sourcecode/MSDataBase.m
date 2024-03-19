classdef MSDataBase
    %MSDATABASE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        DataBaseConnection
        DataBaseFile    (1,1) string
        Features        (:,:) double
        FetchedDataMS1  (:,1) cell
        FetchedDataMS2  (:,1) cell
        Significant     (:,1) logical
        HighFoldChange  (:,1) logical
        FoundInGroup    (:,:) string
        %UI elements
        Window
        Labels
        SelectDBFile
        GenerateDB
        MS1SearchButton
        MS2SearchButton
        mzTolEditfield
        ErrorUnitDropDown
        FilterFragmentationTypeCheckbox
        FilterFragmentationEnergyCheckbox
        FilterSignificantCheckbox
        FilterFoldChangeCheckbox
    end

    methods
        function obj = MSDataBase(CallingApp)
            %MSDATABASE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Features = CallingApp.FeatureData.IdentifierArray;
            IsSignificant = horzcat(CallingApp.FeatureData.SignificantFeature{:});
            IsSignificant = any(IsSignificant <= CallingApp.maxP,2);
            obj.Significant = IsSignificant;
            Fold = max(horzcat(CallingApp.FeatureData.FullFoldChanges{:}),[],2);
            obj.HighFoldChange = Fold >= CallingApp.minFold;
            obj.FoundInGroup = CallingApp.FeatureData.InGroup;
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
            obj.FilterFragmentationTypeCheckbox = uicheckbox(obj.Window,...
                "Value",false,...
                "Position",[60,30,250,20],...
                "Text","Filter fragmentation type",...
                Tooltip="Only consider database MS2 spectra with matching fragmentation type.");
            obj.FilterFragmentationEnergyCheckbox = uicheckbox(obj.Window,...
                "Value",false,...
                "Position",[60,55,250,20],...
                "Text","Filter fragmentation energy",...
                Tooltip="Only consider database MS2 spectra with matching fragmentation energy.");
            obj.FilterSignificantCheckbox = uicheckbox(obj.Window,...
                "Value",false,...
                "Position",[60,80,250,20],...
                "Text","Search only significant features",...
                Tooltip="Only consider significant features in database search.");
            obj.FilterFoldChangeCheckbox = uicheckbox(obj.Window,...
                "Value",false,...
                "Position",[60,105,250,20],...
                "Text","Search only features with high fold change",...
                Tooltip="Only consider features with fold change greater than minimum fold change in database search.");
            drawnow
            obj = obj.CheckStatus;
        end

        function obj=testConnection(obj)          %check for existing SQLlite database
            if isfile(obj.DataBaseFile)
                obj.DataBaseConnection = sqlite(file);
            else
                obj = obj.GenerateDatabase;
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
                obj = obj.LoadDatabase;
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
            obj  = obj.GenerateMassBankDatabase(MassBankRoot);
            close(d)
            obj = obj.CheckStatus;
            uialert(obj.Window,"Local database generation successfull","Database file created","Icon","success")
        end

        function obj = LoadDatabase(obj)
            if isfile(obj.DataBaseFile)
                obj.DataBaseConnection = sqlite(obj.DataBaseFile);
            else
                return
            end
        end

        function DecodedMS2Data = DecodeString(obj,EncodedStrings)
            DecodedMS2Data=cell(size(EncodedStrings,1),1);
            parfor n=1:size(EncodedStrings,1)
                Decoded = matlab.net.base64decode(EncodedStrings(n));
                Decoded = typecast(Decoded,'double');
                DecodedMS2Data{n} = reshape(Decoded,[],3);
            end
        end

        function obj = GenerateMassBankDatabase(obj,rootDir)
            %clear database if one exists
            if ~isempty(obj.DataBaseConnection)
                close(obj.DataBaseConnection)
            end
            % Use the 'dir' function to list all files and folders in the root directory
            fileList = dir(fullfile(rootDir, '**', '*.txt'));
            % Create an empty cell array to store your data
            database = struct('ACCESSION', {}, 'AUTHORS', {}, 'NAME', {}, 'FORMULA', {}, ...
                'CAS_NUMBER', {}, 'EXACT_MASS', {}, 'INSTRUMENT_TYPE', {},'FRAGMENTATION_TYPE', {}, ...
                'PRECURSOR_TYPE', {}, 'PRECURSOR_MZ', {}, 'SPECTRUM', {});
            % Loop through the fileList and read the text files
            parfor i = 1:numel(fileList)
                filename = strcat(fileList(i).folder,'\',fileList(i).name);
                % Read the file into a cell array
                fileLines = importdata(filename, '\n');
                % Initialize variables to store information
                accession = '';
                authors = '';
                name = '';
                formula = '';
                link_cas = '';
                exact_mass = [];
                instrument_type = '';
                fragmentation_type = '';
                focused_ion_type = '';
                focused_ion_mz = [];
                peaks = [];

                % Loop through each line in the file
                for j = 1:numel(fileLines)
                    line = fileLines{j};
                    chNameCaptured = false
                    % Extract the information you need from each line
                    if startsWith(line, 'ACCESSION:')
                        accession = strtrim(strrep(line, 'ACCESSION:', ''));
                    elseif startsWith(line, 'AUTHORS:')
                        authors = strtrim(strrep(line, 'AUTHORS:', ''));
                    elseif ~chNameCaptured && startsWith(line, 'CH$NAME:')
                        name = strtrim(strrep(line, 'CH$NAME:', ''));
                        chNameCaptured = true;  % Set the flag to true to indicate that first CH$NAME has been captured
                    elseif startsWith(line, 'CH$FORMULA:')
                        formula = strtrim(strrep(line, 'CH$FORMULA:', ''));
                    elseif startsWith(line, 'CH$LINK: CAS')
                        link_cas = strtrim(strrep(line, 'CH$LINK: CAS', ''));
                    elseif startsWith(line, 'CH$EXACT_MASS:')
                        exact_mass = str2double(strtrim(strrep(line, 'CH$EXACT_MASS:', '')));
                    elseif startsWith(line, 'AC$INSTRUMENT_TYPE:')
                        instrument_type = strtrim(strrep(line, 'AC$INSTRUMENT_TYPE:', ''));
                    elseif startsWith(line, 'AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE')
                        fragmentation_type = strtrim(strrep(line, 'AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE', ''));
                    elseif startsWith(line, 'MS$FOCUSED_ION: PRECURSOR_TYPE')
                        focused_ion_type = strtrim(strrep(line, 'MS$FOCUSED_ION: PRECURSOR_TYPE', ''));
                    elseif startsWith(line, 'MS$FOCUSED_ION: PRECURSOR_M/Z')
                        focused_ion_mz = str2double(strtrim(strrep(line, 'MS$FOCUSED_ION: PRECURSOR_M/Z', '')));
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
                                peaks = [peaks; str2double(parts(1)), str2double(parts(2)), str2double(parts(3))];
                            end
                        end
                    end
                end
                %check inputs
                if isempty(focused_ion_mz)
                    focused_ion_mz = NaN;
                end
                if isempty(exact_mass)
                    exact_mass = NaN;
                end
                if isempty(peaks)
                    peaks = [0 0 0];
                end
                %decode peaks to string
                peaks=reshape(peaks,1,[]);
                peaks=typecast(peaks,'uint8');
                peaks=matlab.net.base64encode(peaks);
                % Add the extracted data to the database structure
                entry = struct('ACCESSION', accession, 'AUTHORS', authors, 'NAME', name, ...
                    'FORMULA', formula, 'CAS_NUMBER', link_cas, 'EXACT_MASS', exact_mass, ...
                    'INSTRUMENT_TYPE', instrument_type, 'FRAGMENTATION_TYPE', fragmentation_type, 'PRECURSOR_TYPE', focused_ion_type, ...
                    'PRECURSOR_MZ', focused_ion_mz, 'SPECTRUM', peaks);

                database = [database, entry];
            end
            database=struct2table(database);
            obj.DataBaseConnection = sqlite(obj.DataBaseFile,"create");
            sqlwrite(obj.DataBaseConnection,"MassBankData",database);
        end

        function obj = MassBankMS1Query(obj)
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
            
            query = ['SELECT DISTINCT NAME, ' ...
                    '	FORMULA, ' ...
                    '	CAS_NUMBER, ' ...
                    '	EXACT_MASS, ' ...
                    'FROM MassBankData ' ...
                    'WHERE EXACT_MASS <= '];
            query = append(query,convertStringsToChars(MZmax + " AND EXACT_MASS >= " + MZmin));
            databasefile = obj.DataBaseFile;
            QueryResults = cell(size(query));
            parfor n=1:length(QueryMasses)
                local_connection = sqlite(databasefile)
                result = fetch(local_connection, query{n});
                %calculate difference in ppm
                diff = ((abs(QueryMasses(n)-result.EXACT_MASS))./QueryMasses(n))*10^6;
                result.DELTA_ppm = round(diff,2);
                result = sortrows(result,"DELTA_ppm","ascend");
                close(local_connection);
                if ~isempty(result)
                    %copy query to all features with same mass
                    QueryResults{n}=result;
                end
            end
            %sort data to feature
            out = cell(length(obj.Features(:,1)),1);
            for n=1:size(QueryResults,1)
                idx = obj.Features(:,1) == QueryMasses(n);
                out(idx) = QueryResults(n);
            end
            out(~Modifier) = {};
            obj.FetchedDataMS1 = out; 
        end
        
        function obj = MassBankMS2Query(obj)
            %Query time
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
            
            query = ['SELECT DISTINCT NAME, ' ...
                    '	FORMULA, ' ...
                    '	CAS_NUMBER, ' ...
                    '	EXACT_MASS, ' ...
                    '	PRECURSOR_TYPE ' ...
                    '	FRAGMENTATION_TYPE ' ...
                    '	SPECTRUM ' ...
                    'FROM MassBankData ' ...
                    'WHERE EXACT_MASS <= '];
            query = append(query,convertStringsToChars(MZmax + " AND EXACT_MASS >= " + MZmin));
            databasefile = obj.DataBaseFile;
            QueryResults = cell(size(query));

            parfor n=1:length(QueryMasses)
                local_connection = sqlite(databasefile)
                result = fetch(local_connection, query{n});
                if ~isempty(result)
                    %copy query to all features with same mass
                    QueryResults{n}=result;
                end
                close(local_connection);
            end
            %sort data to feature
            out = cell(length(obj.Features(:,1)),1);
            for n=1:size(QueryResults,1)
                idx = obj.Features(:,1) == QueryMasses(n);
                out(idx) = QueryResults(n);
            end
            out(~Modifier) = {};
            % calculate spectral match score

            obj.FetchedDataMS2 = out; 
        end

        function obj = CheckStatus(obj)
            if isempty(obj.DataBaseFile)
                obj.MS1SearchButton.Enable = "off";
                obj.MS2SearchButton.Enable = "off";
                obj.mzTolEditfield.Enable = "off";
                obj.ErrorUnitDropDown.Enable = "off";
                obj.FilterFragmentationTypeCheckbox.Enable = "off";
                obj.FilterFragmentationEnergyCheckbox.Enable = "off";
                obj.FilterSignificantCheckbox.Enable = "off";
                obj.FilterFoldChangeCheckbox.Enable = "off";
            else
                obj.MS1SearchButton.Enable = "on";
                obj.MS2SearchButton.Enable = "on";
                obj.mzTolEditfield.Enable = "on";
                obj.ErrorUnitDropDown.Enable = "on";
                obj.FilterFragmentationTypeCheckbox.Enable = "on";
                obj.FilterFragmentationEnergyCheckbox.Enable = "on";
                obj.FilterSignificantCheckbox.Enable = "on";
                obj.FilterFoldChangeCheckbox.Enable = "on";
            end

            if isempty(obj.DataBaseConnection)
                obj.MS1SearchButton.Enable = "off";
                obj.MS2SearchButton.Enable = "off";
                obj.mzTolEditfield.Enable = "off";
                obj.ErrorUnitDropDown.Enable = "off";
                obj.FilterFragmentationTypeCheckbox.Enable = "off";
                obj.FilterFragmentationEnergyCheckbox.Enable = "off";
                obj.FilterSignificantCheckbox.Enable = "off";
                obj.FilterFoldChangeCheckbox.Enable = "off";
            else
                obj.MS1SearchButton.Enable = "on";
                obj.MS2SearchButton.Enable = "on";
                obj.mzTolEditfield.Enable = "on";
                obj.ErrorUnitDropDown.Enable = "on";
                obj.FilterFragmentationTypeCheckbox.Enable = "on";
                obj.FilterFragmentationEnergyCheckbox.Enable = "on";
                obj.FilterSignificantCheckbox.Enable = "on";
                obj.FilterFoldChangeCheckbox.Enable = "on";
            end
            drawnow
        end


    end
end

