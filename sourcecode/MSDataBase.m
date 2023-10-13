classdef MSDataBase
    %MSDATABASE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        DataBaseConnection
        DataBaseFile (1,1) string
        Features (:,:) double
        FetchedData (1,:) cell
        Significance (1,:) cell
        FoundInGroup (:,:) string
        mzTol       (1,1) double {mustBeFinite} = 0.001
        ErrorUnit   (1,1) string {mustBeMember(ErrorUnit,["Da","ppm"])} = "Da"
        FragmentationType (1,1) logical = false
        %UI elements
        Window
        Labels
        MS1SearchButton
        MS2SearchButton
        RebuildMassBankDataBaseButton
        mzTolEditfield
        ErrorUnitDropDown
        FilterFragmentationTypeCheckbox
    end

    methods
        function obj = MSDataBase(CallingApp)
            %MSDATABASE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Features = CallingApp.FeatureData.IdentifierArray;
            obj.Significance = CallingApp.FeatureData.SignificantFeature;
            obj.FoundInGroup = CallingApp.FeatureData.InGroup;

            %build UIelements
            obj.Window = uifigure("WindowStyle","normal",...
                "Name","Database Searcher",...
                "Position",[100,100,700,480],...
                "NumberTitle","off");
            obj.MS1SearchButton = uibutton(obj.Window,"push",...
                "Position",[535,420,155,40],...
                Text="MS1 search");
            obj.mzTolEditfield = uieditfield(obj.Window,"numeric", ...
                "Position",[435,440,40,20], ...
                "Limits",[0,Inf], ...
                "LowerLimitInclusive","off",...
                "Value",obj.mzTol,...
                Tooltip="Lowest mass tolerance to identify features between different groups as the same feature.");
            obj.Labels(2) = uilabel(obj.Window,"Text","Allowed mass tolerance","Position",[280,440,150,20],HorizontalAlignment="right");

            obj.ErrorUnitDropDown = uidropdown(obj.Window, ...
                "Items",["ppm","Da"],...
                ItemsData=["ppm","Da"], ...
                Value=obj.ErrorUnit,...
                Position=[475,440,50,20],...
                Tooltip="Specify whether a fixed mz error in Da or a relative error in ppm is used.");
            drawnow
            %check for existing SQLlite database
            if isfile(obj.DataBaseFile)
                obj.DataBaseConnection = sqlite(file);
            else
                obj.GenerateDatabase
            end
        end

        function obj = GenerateDatabase(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [file,path] = uiputfile('*.db','Select Save Location and Filename for database file',"DataBase.db");
            obj.DataBaseFile=fullfile(path,file);
            MassBankRoot = uigetdir(pwd,'Open MassBank data main folder');
            d =  uiprogressdlg(obj.Window,'Title','Building Database',...
            'Indeterminate','on');
            drawnow
            obj  = obj.GenerateMassBankDatabase(MassBankRoot);
            close(d)
        end

        function outputArg = LoadDatabase(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            files = matlab.apputil.getInstalledAppInfo;
            id = [files.id] == "AriumMSAPP";
            path = files(id).location;
            file= fullfile(path,"SpectralDataBase.db");
            if isfile(file)
                obj.DataBaseConnection = sqlite(file);
            else
                return
            end
        end

        function DecodedMS2Data = DecodeString(obj,EncodedStrings)
                DecodedMS2Data=cell(size(EncodedStrings,1),1);
            parfor n=size(EncodedStrings,1)
                Decoded = matlab.net.base64decode(EncodedStrings(n));
                Decoded = typecast(Decoded,'double');
                DecodedMS2Data{n} = reshape(Decoded,[],3);
            end
        end

        function obj = GenerateMassBankDatabase(obj,rootDir)
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

        function obj=MassBankQuery(obj,mzTolerance,ErrorUnit,FragmentationType)

            query = ['SELECT NAME, ' ...
                '	FORMULA, ' ...
                '	CAS_NUMBER, ' ...
                '	EXACT_MASS, ' ...
                '	PRECURSOR_TYPE, ' ...
                '	SPECTRUM ' ...
                'FROM SpectralDataBase ' ...
                'WHERE FRAGMENTATION_TYPE IN (FragType) ' ...
                '	AND EXACT_MASS - QueryMZ <= tolerance'];

            %% Execute query and fetch results
            obj.FetchedData = fetch(conn,query);
        end

    end
end

