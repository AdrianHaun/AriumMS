function [instrument,ionization,energy,Spectrum] = readHMDBSpectraFile(SpectraFile)

Spectraldoc = xmlread(SpectraFile);

instrument = Spectraldoc.getElementsByTagName('instrument-type');
instrument = string(instrument.item(0).getTextContent);

ionization = Spectraldoc.getElementsByTagName('ionization-mode');
ionization = string(ionization.item(0).getTextContent);

energy = Spectraldoc.getElementsByTagName('collision-energy-voltage');
energy = str2double(energy.item(0).getTextContent);

Spectra = Spectraldoc.getElementsByTagName('ms-ms-peak');

Spectrum = zeros(Spectra.getLength,2);

for n = 0:Spectra.getLength-1
    Peak = Spectra.item(n);
    mz = Peak.getElementsByTagName('mass-charge');
    Spectrum(n+1,1) = str2double(mz.item(0).getTextContent);

    Int = Peak.getElementsByTagName('intensity');
    Spectrum(n+1,2) = str2double(Int.item(0).getTextContent);
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
if isempty(Spectrum)
    Spectrum = [0 0];
end

Spectrum = reshape(Spectrum,1,[]);
Spectrum=typecast(Spectrum,'uint8');
Spectrum=matlab.net.base64encode(Spectrum);

end