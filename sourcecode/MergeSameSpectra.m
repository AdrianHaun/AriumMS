function MergedSpectra = MergeSameSpectra(Spectra,Scores)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%return if only one spectra present
if size(Scores,1)<=2
    MergedSpectra = Spectra;
    return
end

Scores(Scores(:,1)<850,:)=[];
SpectraIndex = Scores(:,2:3);

if isempty(Scores)
    MergedSpectra = Spectra;
    return
end

nSpectra = unique(SpectraIndex);
SortedSpectraIdx = cell(size(nSpectra));
for n=1:length(nSpectra)
    %initial index
    Query = nSpectra(n);
    id = any(SpectraIndex==Query,2);
    Cluster = unique(SpectraIndex(id,:));
    HasNewSpectra=true;
    while HasNewSpectra
        id = any(ismember(SpectraIndex,Cluster),2);
        newSpec = unique(SpectraIndex(id,:));
        if any(~ismember(newSpec,Cluster))
            Cluster = unique(SpectraIndex(id,:));
            HasNewSpectra = true;
        else
            HasNewSpectra = false;
        end
    end
    if size(Cluster,1) > size(Cluster,2)
        SortedSpectraIdx{n} = Cluster;
    else
        SortedSpectraIdx{n} = Cluster';
    end
end

% Initialize an empty cell array to store unique content
uniqueCluster = {};

% Iterate through the cell array
for i = 1:numel(SortedSpectraIdx)
    % Check if the current cell content exists in uniqueCluster
    isDuplicate = false;
    for j = 1:numel(uniqueCluster)
        if isequal(SortedSpectraIdx{i}, uniqueCluster{j})
            isDuplicate = true;
            break;
        end
    end
    
    % If the content is not a duplicate, add it to uniqueCluster
    if ~isDuplicate
        uniqueCluster{end+1,1} = SortedSpectraIdx{i};
    end
end

%prepare Spectra Matrix
Masses = Spectra(:,1);
Spectra(:,1) = [];
MergedSpectra = cell(1,numel(uniqueCluster));
for n=1:numel(uniqueCluster)
    MergedSpectra{n} = mean(Spectra(:,uniqueCluster{n}),2);
end

MergedSpectra = [Masses, horzcat(MergedSpectra{:})];