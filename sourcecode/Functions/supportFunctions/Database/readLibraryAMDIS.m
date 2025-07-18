function entryStringCell = readLibraryAMDIS(fileLocation)

arguments
    fileLocation (1,1) string {mustBeFile}
end

%% import TXT file
dataLines = [1, Inf];

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ":";

% Specify column names and types
opts.VariableNames = ["NAME", "MinAlanineMCFDerivativeMainPeak", "Var3"];
opts.SelectedVariableNames = ["NAME", "MinAlanineMCFDerivativeMainPeak"];
opts.VariableTypes = ["string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["NAME", "MinAlanineMCFDerivativeMainPeak", "Var3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["NAME", "MinAlanineMCFDerivativeMainPeak", "Var3"], "EmptyFieldRule", "auto");

% Import the data
fileStringArray = readmatrix(fileLocation, opts);
%% process fileString

%split each entry
splitIndex = find(ismissing(fileStringArray(:,1)));
splitIndex = diff([1;splitIndex]);
splitIndex(1) = splitIndex(1) + 1;

entryStringCell = mat2cell(fileStringArray,splitIndex);


parfor iEntry = 1:height(entryStringCell)
    currentEntry = entryStringCell{iEntry,1};
    %remove RT from Name
    currentEntry(1,2) = extractAfter(currentEntry(1,2),wildcardPattern+"min ");
      
    index = find(currentEntry(:,1)=="NUM PEAKS")+1; %find spectrum start index
    spectrumString = currentEntry(index:end-1,1);   %isolate spectrum data
    currentEntry([4:5,7:end-2],:) = []; %remove unnecessary fields

    %process spectrum data
    spectrumString = strip(spectrumString,"both"); %remove whitespace
    spectrumString = join(spectrumString'); 
    spectrumString = split(spectrumString,") (");
    %remove remaining brackets 
    spectrumString = strip(spectrumString,"(");
    spectrumString = strip(spectrumString,")");
    spectrumString = replace(spectrumString,whitespacePattern," ");
    spectrumString = strip(spectrumString,"both"); %remove possible leading or terminating whitespace
    spectrumString = split(spectrumString," ",2);
    %rescale intensities
    spectrumString = str2double(spectrumString);
    spectrumString(:,2) = spectrumString(:,2)./max(spectrumString(:,2));
    %store max intensity fragment as RefIon
    referenceMass = spectrumString(spectrumString(:,2)==1,1);

    currentEntry(5,:) = ["Ref.Mass", string(referenceMass(1))];
    %store as single string
    spectrumString = reshape(spectrumString,1,[]);
    spectrumString = join(string(spectrumString),"-");
    currentEntry(end,:) = ["SPECTRUM", spectrumString];
    %store processed entry
    entryStringCell{iEntry,1} = currentEntry;
end
