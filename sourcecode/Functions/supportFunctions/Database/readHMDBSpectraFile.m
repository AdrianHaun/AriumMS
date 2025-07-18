function [instrument,ionization,energy,spectrum] = readHMDBSpectraFile(spectraFile)
%% readHMDBSpectraFile extracts relevant MS and MSMS information from HMDB spectra file (.XML)

arguments
    spectraFile (1,1) string {mustBeFile}
end

spectraldoc = xmlread(spectraFile);

instrument = spectraldoc.getElementsByTagName('instrument-type');
instrument = string(instrument.item(0).getTextContent);

ionization = spectraldoc.getElementsByTagName('ionization-mode');
ionization = string(ionization.item(0).getTextContent);

energy = spectraldoc.getElementsByTagName('collision-energy-voltage');
energy = str2double(energy.item(0).getTextContent);

spectrumArray = spectraldoc.getElementsByTagName('ms-ms-peak');

spectrum = zeros(spectrumArray.getLength,2);

for iSpectrum = 0:spectrumArray.getLength-1
    peak = spectrumArray.item(iSpectrum);
    mass = peak.getElementsByTagName('mass-charge');
    spectrum(iSpectrum+1,1) = str2double(mass.item(0).getTextContent);

    intensity = peak.getElementsByTagName('intensity');
    spectrum(iSpectrum+1,2) = str2double(intensity.item(0).getTextContent);
end

if isempty(instrument) || strcmp(instrument,"")
    instrument = "Unknown";
end
if isempty(ionization) || strcmp(ionization,"")
    ionization = "Unknown";
end
if isempty(energy) || isnan(energy)
    energy = 0;
end
if isempty(spectrum)
    spectrum = [0 0];
end

spectrum = reshape(spectrum,1,[]);
spectrum = typecast(spectrum,'uint8');
spectrum = matlab.net.base64encode(spectrum);