function CompundScores = InnerFeatScores(MSnSpectraMatrix)
Intensities = MSnSpectraMatrix(:,2:end);
%mass weighted intensities
Intensities = Intensities.*MSnSpectraMatrix(:,1);

%scale each spectrum
Intensities = Intensities./max(Intensities);

NumSpectra = size(Intensities,2);
if NumSpectra == 1
    CompundScores = [999,1,1];
else
SpectraIndex = 1:1:NumSpectra;
Combis = nchoosek(SpectraIndex,2);
CompundScores = zeros(size(Combis,1),1);
for n = 1:size(Combis,1)
    Spectra1 = Intensities(:,Combis(n,1))';
    Spectra2 = Intensities(:,Combis(n,2))';
    idx = Spectra1 == 0 & Spectra2 == 0;
    Spectra1(idx) = [];
    Spectra2(idx) = [];
    CompundScores(n,1) = CompositScore(Spectra1,Spectra2);
end
CompundScores = [CompundScores,Combis];
end