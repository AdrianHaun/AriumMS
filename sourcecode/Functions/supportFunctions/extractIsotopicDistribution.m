function [isoDistribution, chargeState, elementHits] = extractIsotopicDistribution(targetMass, spectrum, varargin)
%EXTRACTISOTOPICDISTRIBUTION Extracts isotopic distribution and estimates charge state.
%
%   [isoMasses, isoIntensities, chargeState] = extractIsotopicDistribution(targetMass, spectrum)
%
% Optional parameters:
%   'Window'    - scalar, range above targetMass (default: 5)
%   'Tolerance' - scalar, tolerance in Da for grouping (default: 0.01)
%
% Outputs:
%   isoDistribution - masses of detected isotopic peaks (column1), and normalized
%                       intensities (column2)
%   chargeState    - estimated integer charge state

    % Defaults
    window = 5;
    tolerance = 0.01;

    % Parse optional arguments
    for i = 1:2:length(varargin)
        param = varargin{i};
        value = varargin{i+1};
        switch lower(param)
            case 'window'
                window = value;
            case 'tolerance'
                tolerance = value;
            otherwise
                error('Unknown parameter: %s', param);
        end
    end

    % Adjusted mass window
    minMass = targetMass - 0.1;
    maxMass = targetMass + window + 0.1;

    % Subset the spectrum
    mask = spectrum(:,1) >= minMass & spectrum(:,1) <= maxMass;
    subSpec = spectrum(mask, :);

    if isempty(subSpec)
        isoDistribution = [];
        chargeState = NaN;
        return;
    end

    % Sort by mass
    [~, sortIdx] = sort(subSpec(:,1));
    subSpec = subSpec(sortIdx, :);
    masses = subSpec(:,1);
    intensities = subSpec(:,2);

    % Initialize cluster detection
    clusters = {};
    used = false(size(masses));

    % Loop over all peaks to build clusters
    for i = 1:length(masses)
        if used(i)
            continue;
        end
        thisMass = masses(i);
        clusterMasses = thisMass;
        clusterIntensities = intensities(i);
        used(i) = true;

        % Try to find subsequent isotopes
        nextMass = thisMass + 1;
        while true
            % Find closest peak within tolerance
            delta = abs(masses - nextMass);
            [minDelta, idx] = min(delta);
            if minDelta <= tolerance && ~used(idx)
                clusterMasses(end+1) = masses(idx); %#ok<AGROW>
                clusterIntensities(end+1) = intensities(idx); %#ok<AGROW>
                used(idx) = true;
                nextMass = nextMass + 1;
            else
                break;
            end
        end

        % Store this cluster
        clusters{end+1} = struct('masses', clusterMasses, 'intensities', clusterIntensities); %#ok<AGROW>
    end

    % Choose cluster whose first mass is closest to targetMass
    minDist = inf;
    bestIdx = 0;
    for k = 1:numel(clusters)
        dist = abs(clusters{k}.masses(1) - targetMass);
        if dist < minDist
            minDist = dist;
            bestIdx = k;
        end
    end

    if bestIdx == 0
        isoDistribution = [];
        chargeState = NaN;
        return;
    end

    % Get selected cluster
    selectedMasses = clusters{bestIdx}.masses(:);
    selectedIntensities = clusters{bestIdx}.intensities(:);

    % Estimate charge state by spacing
    if numel(selectedMasses) >= 2
        spacings = diff(selectedMasses);
        meanSpacing = mean(spacings);
        chargeEstimate = round(1 / meanSpacing);
        % Sanity check
        if chargeEstimate < 1
            chargeEstimate = 1;
        end
    else
        chargeEstimate = 1; % default
    end

    % Normalize intensities
    if max(selectedIntensities) > 0
        normIntensities = selectedIntensities / max(selectedIntensities);
    else
        normIntensities = selectedIntensities;
    end
    
    % Output
    isoDistribution = [selectedMasses,normIntensities];
    chargeState = chargeEstimate;
    elementHits = detectIsotopicElements(selectedMasses, normIntensities, chargeState);
end
