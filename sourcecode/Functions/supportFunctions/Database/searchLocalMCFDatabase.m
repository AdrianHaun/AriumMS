function results = searchLocalMCFDatabase(FeatureStruct,dataBaseTable)

compoundStringArray = strings(FeatureStruct.dataSize,1);
formulaStringArray = strings(FeatureStruct.dataSize,1);
CASStringArray = strings(FeatureStruct.dataSize,1);

scoreArray =  zeros(FeatureStruct.dataSize,1);


%decode database spectra
dataBaseSpectra = decodeSpectra(dataBaseTable.SPECTRUM);
%% search database

parfor iFeature = 1:FeatureStruct.dataSize

    featureSpectrum = FeatureStruct.feature(iFeature).spectrumMS2;
    featureSpectrum(:,1) = round(featureSpectrum(:,1));

    allScores = zeros(height(dataBaseTable),1);

    for jEntry = 1:height(dataBaseSpectra)
        databaseSpectrum = dataBaseSpectra(jEntry);
        spectraCells = [featureSpectrum;databaseSpectrum];
        %align spectra
        alignedSpectra = alignSpectra(spectraCells,"normal","low","true",0.05);
        %composite score
        compositScore = calculateCompositeScore(alignedSpectra(:,2),alignedSpectra(:,3));
        allScores(jEntry) = compositScore;
    end
    [score,ID] = max(allScores);
    scoreArray(iFeature) = score;
    compoundStringArray(iFeature) = dataBaseTable.NAME(ID);
    formulaStringArray(iFeature) = dataBaseTable.FORMULA(ID);
    CASStringArray(iFeature) = dataBaseTable.CAS_NUMBER(ID);
end
results = table(compoundStringArray,formulaStringArray,CASStringArray,scoreArray,'VariableNames',["Name","Formula","CAS","Score"]);
