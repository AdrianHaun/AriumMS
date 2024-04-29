function result = ResampleMS2Spectra(Spectrum)

% Tolerance value for grouping masses
tolerance = 0.05;

% Sort the data based on mz values
sortedData = sortrows(Spectrum, 1);

% Initialize variables to store grouped masses and intensities
groupedMasses = [];
groupedIntensities = [];
currentGroup = [];
currentSum = [];
count = 0;

% Iterate through the sorted data
for i = 1:size(sortedData, 1)
    mz = sortedData(i, 1);
    intensity = sortedData(i, 2);
    % If the current mz is within tolerance of the previous mz, add it to the group
    if isempty(currentGroup) || abs(mz - currentGroup(end)) <= tolerance
        currentGroup = [currentGroup, mz];
        currentSum = [currentSum, intensity];
        count = count + 1;
    else
        % Take mz with highest intensity in the group
        [~,idx] = max(currentSum,[],"all");
        averageMass = currentGroup(idx);
        currentSum = sum(currentSum);
        % Store the mass and corresponding summed intensity
        groupedMasses = [groupedMasses; averageMass];
        groupedIntensities = [groupedIntensities; currentSum];
        
        % Reset variables for the next group
        currentGroup = mz;
        currentSum = intensity;
        count = 1;
    end
end

% Handle the last group
[~,idx] = max(currentSum,[],"all");
averageMass = currentGroup(idx);
currentSum = sum(currentSum);
% Store the mass and corresponding summed intensity
groupedMasses = [groupedMasses; averageMass];
groupedIntensities = [groupedIntensities; currentSum];

%Combine the result
result = [groupedMasses, groupedIntensities];
end