function a = augmentFeatures(a, b)

%preallocate id storage
id = cell(length(b),1);

tolerance = 0.05; % Define the tolerance

% Loop through each entry in struct b
parfor i = 1:length(b)
    % Find the index of mass_measured in a that is within the specified tolerance
    id{i} = find(abs([a.mass_measured] - b(i).mass_measured) <= tolerance);
end

% Initialize two new structs
% Get fieldnames of struct b
fieldNames = fieldnames(b);

% Initialize two new structs with the same fieldnames
b_no_match = struct();
b_with_match = struct();

% Preallocate the structs with the same fieldnames
for k = 1:length(fieldNames)
    b_no_match.(fieldNames{k}) = [];
    b_with_match.(fieldNames{k}) = [];
end


% Loop through each entry in struct b to split into two structs
for i = 1:length(b)
    if isempty(id{i})
        b_no_match(end+1) = b(i); % Add to b_empty if id is empty
    else
        b_with_match(end+1) = b(i); % Add to b_non_empty if id has values
    end
end
% Remove the first empty entry from the structs
b_no_match(1) = [];
b_with_match(1) = [];

id = id(~cellfun('isempty', id)); % Remove empty cells in id

% add b_with_match to a based on spectrumMS2
for i = 1:height(id)
    id_of_a = id{i};
    check = numel(id_of_a);

    spectrum_of_b = {b(i).spectrumMS2};

    switch check
        case 1 % case 1 only one entry
            spectrum_of_a = {a(id_of_a).spectrumMS2};
            alignedSpectra = alignSpectra([spectrum_of_b;spectrum_of_a],"normal","low","false");
            %calculate score
            score = calculateCompositeScore(alignedSpectra(:,2),alignedSpectra(:,3));

            if score > 700
                % append all entries of b inside entry of a
                
            else
                continue
            end

        otherwise % case 2: multiple entries

    end


end
end