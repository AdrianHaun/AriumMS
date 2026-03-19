function data = readMSfile_pwiz(filename,dllPath)
%READMZML_PWIZ Read mzML using ProteoWizard (pwiz) .NET bindings into a MATLAB struct.

if nargin < 1 || isempty(filename)
    error('You must provide an MS filename.');
end
if ~isfile(filename)
    error('File not found: %s', filename);
end

if nargin < 2
    dllPath = '';
end

% Load pwiz .NET assembly
if ~isempty(dllPath)
    if ~isfile(dllPath)
        error('pwiz .NET DLL not found: %s', dllPath);
    end
    NET.addAssembly(dllPath);
end

% Import namespaces
try
    import pwiz.CLI.msdata.*;
    import pwiz.CLI.cv.*;
catch
    error(['Unable to find pwiz.CLI.msdata or pwiz.CLI.cv namespaces. ' ...
        'Make sure pwiz_bindings_cli.dll is on the path or provided as pwizDllPath.']);
end

% Open the ms file
msd = MSDataFile(filename);
spectrumList  = msd.run.spectrumList;
nScans   = spectrumList.size();

spectra(nScans) = struct( ...
    'msLevel',           [], ...
    'rt',                [], ...
    'rawScan',           [], ...
    'polarity',          '', ...
    'precursorMz',       [], ...
    'precursorCharge',   [], ...
    'precursorMass',     [], ...
    'collisionEnergy',   [], ...
    'fragmentationType', '', ...
    'basePeakMz',        [], ...
    'basePeakIntensity', [], ...
    'totalIonCurrent',   [], ...
    'dataType',          '', ...
    'processedScan',     [], ...
    'centroidedScan',    []);

for iScan = 0:nScans-1
    spectrum = spectrumList.spectrum(iScan, true);   % with binary data

    % --- m/z and intensity arrays ---
    mzArray  = spectrum.getMZArray();
    intArray = spectrum.getIntensityArray();

    if isempty(mzArray)
        mz = [];
    else
        mz = double(mzArray.data.Storage);
    end

    if isempty(intArray)
        intensity = [];
    else
        intensity = double(intArray.data.Storage);
    end

    scanData = [mz;intensity]';


    % --- MS level ---
    msLevel = NaN;
    try
        msLevelParam = spectrum.cvParam(CVID.MS_ms_level);
        if ~msLevelParam.empty()
            rawVal = msLevelParam.value;
            try
                v = str2double(char(rawVal.ToString()));
                if ~isnan(v), msLevel = v; end
            catch
            end
            if isnan(msLevel)
                try
                    vNet = System.Convert.ToDouble(rawVal);
                    msLevel = double(vNet);
                catch
                end
            end
        end
    catch
    end

     % --- Data type: profile or centroid ---
        dataType = '';
        try
            % Profile spectrum?
            try
                profParam = spectrum.cvParam(CVID.MS_profile_spectrum);
            catch
                profParam = [];
            end
            % Centroided spectrum?
            try
                centParam = spectrum.cvParam(CVID.MS_centroid_spectrum);
            catch
                centParam = [];
            end

            if ~isempty(profParam) && ~profParam.empty()
                dataType = 'profile';
            elseif ~isempty(centParam) && ~centParam.empty()
                dataType = 'centroid';
            end
        catch
        end

    % --- RT (seconds) ---
    rt = NaN;
    scans = spectrum.scanList.scans;
    if scans.Count > 0
        scan0 = scans.Item(int32(0));
        try
            rtParam = scan0.cvParam(CVID.MS_scan_start_time);
            if ~rtParam.empty()
                rt = double(rtParam.timeInSeconds());
            end
        catch
        end
    end

    % --- Polarity ---
    polarity = '';
    try
        posParam = spectrum.cvParam(CVID.MS_positive_scan);
        negParam = spectrum.cvParam(CVID.MS_negative_scan);
        if exist('posParam','var') && ~posParam.empty()
            polarity = '+';
        elseif exist('negParam','var') && ~negParam.empty()
            polarity = '-';
        end
    catch
    end

    % --- Precursor info ---
    precursorMz       = NaN;
    precursorCharge   = NaN;
    collisionEnergy   = NaN;
    fragmentationType = '';

    try
        precursors = spectrum.precursors;
        if precursors.Count > 0
            prec = precursors.Item(int32(0));

            % Selected ion
            selIons = prec.selectedIons;
            if selIons.Count > 0
                sel = selIons.Item(int32(0));

                % precursor m/z
                try pMzParam = sel.cvParam(CVID.MS_selected_ion_m_z); catch, pMzParam = []; end
                if ~isempty(pMzParam) && ~pMzParam.empty()
                    rawVal = pMzParam.value;
                    try
                        v = str2double(char(rawVal.ToString()));
                        if ~isnan(v), precursorMz = v; end
                    catch
                    end
                    if isnan(precursorMz)
                        try
                            vNet = System.Convert.ToDouble(rawVal);
                            precursorMz = double(vNet);
                        catch
                        end
                    end
                end

                % charge
                try zParam = sel.cvParam(CVID.MS_charge_state); catch, zParam = []; end
                if ~isempty(zParam) && ~zParam.empty()
                    rawVal = zParam.value;
                    try
                        v = str2double(char(rawVal.ToString()));
                        if ~isnan(v), precursorCharge = v; end
                    catch
                    end
                    if isnan(precursorCharge)
                        try
                            vNet = System.Convert.ToDouble(rawVal);
                            precursorCharge = double(vNet);
                        catch
                        end
                    end
                end
            end

            % Activation
            act = prec.activation;

            % collision energy
            try
                ceParam = act.cvParam(CVID.MS_collision_energy);
                if ~ceParam.empty()
                    rawVal = ceParam.value;
                    try
                        v = str2double(char(rawVal.ToString()));
                        if ~isnan(v), collisionEnergy = v; end
                    catch
                    end
                    if isnan(collisionEnergy)
                        try
                            vNet = System.Convert.ToDouble(rawVal);
                            collisionEnergy = double(vNet);
                        catch
                        end
                    end
                end
            catch
            end

            % fragmentation type
            fragType = '';
            fragCVIDs = { ...
                'MS_collision_induced_dissociation', ...
                'MS_higher_energy_collisional_dissociation', ...
                'MS_electron_transfer_dissociation', ...
                'MS_electron_capture_dissociation'};

            for f = 1:numel(fragCVIDs)
                try
                    cvField = CVID.(fragCVIDs{f});
                    p = act.cvParam(cvField);
                    if ~p.empty()
                        fragType = char(p.name);
                        break;
                    end
                catch
                end
            end
            fragmentationType = fragType;
        end
    catch
    end


    % --- Base peak and TIC ---
    basePeakMz        = NaN;
    basePeakIntensity = NaN;
    totalIonCurrent   = NaN;

    try
        bpMzParam = spectrum.cvParam(CVID.MS_base_peak_m_z);
    catch
        bpMzParam = [];
    end
    if ~isempty(bpMzParam) && ~bpMzParam.empty()
        rawVal = bpMzParam.value;
        try
            v = str2double(char(rawVal.ToString()));
            if ~isnan(v), basePeakMz = v; end
        catch
        end
    end

    try
        bpIntParam = spectrum.cvParam(CVID.MS_base_peak_intensity);
    catch
        bpIntParam = [];
    end
    if ~isempty(bpIntParam) && ~bpIntParam.empty()
        rawVal = bpIntParam.value;
        try
            v = str2double(char(rawVal.ToString()));
            if ~isnan(v), basePeakIntensity = v; end
        catch
        end
    end

    try
        ticParam = spectrum.cvParam(CVID.MS_total_ion_current);
    catch
        ticParam = [];
    end
    if ~isempty(ticParam) && ~ticParam.empty()
        rawVal = ticParam.value;
        try
            v = str2double(char(rawVal.ToString()));
            if ~isnan(v), totalIonCurrent = v; end
        catch
        end
    end

    spectra(iScan+1).msLevel           = msLevel;
    spectra(iScan+1).rt                = rt;
    spectra(iScan+1).rawScan           = scanData;
    spectra(iScan+1).polarity          = polarity;
    spectra(iScan+1).precursorMz       = precursorMz;
    spectra(iScan+1).precursorCharge   = precursorCharge;
    spectra(iScan+1).collisionEnergy   = collisionEnergy;
    spectra(iScan+1).fragmentationType = fragmentationType;
    spectra(iScan+1).basePeakMz        = basePeakMz;
    spectra(iScan+1).basePeakIntensity = basePeakIntensity;
    spectra(iScan+1).totalIonCurrent   = totalIonCurrent;
    spectra(iScan+1).dataType          = dataType;
end

% Run-level start timestamp
startTimeStamp = '';
try
    if ~isempty(msd.run.startTimeStamp)
        startTimeStamp = char(msd.run.startTimeStamp.ToString());
    end
catch
end

% remove empty scans in spectra struct
    % Determine which spectra have non-empty rawScan, msLevel and RT
    isEmptyScan = true(1,numel(spectra));
    for k = 1:numel(spectra)
        s = spectra(k);
        % Consider a scan empty if rawScan is empty or missing, RT is NaN or msLevel is NaN
        hasRaw = ~isempty(s.rawScan);
        hasRT = ~isnan(s.rt);
        hasMsL = ~isnan(s.msLevel);
        if hasRaw && hasRT && hasMsL
            isEmptyScan(k) = false;
        end
    end
    spectra = spectra(~isEmptyScan);

%update number of scans
nScans = numel(spectra);

data = struct();
data.file           = filename;
data.startTimeStamp = startTimeStamp;
data.nSpectra       = nScans;
data.spectra        = spectra;
end