function Report = evaluateDataBaseSearch(Report,results,FeatureStruct)

results = results(:,[1,4]);

timesReport = Report.RT;
timesFeatures = [FeatureStruct.retentionTime]';


foundTime = zeros(height(Report),1);
foundFeatID = foundTime;
foundScore = foundTime;
foundCompund = strings(height(Report),1);
%search each Report entry for matching Feature
for iEntry = 1:height(Report)
    [~,minID] = min(abs(timesFeatures - timesReport(iEntry)));
    foundTime(iEntry) = timesFeatures(minID);
    foundFeatID(iEntry) = minID;
    foundScore(iEntry) = results.Score(minID);
    foundCompund(iEntry) = results.Name(minID);
end
%store results in Report
Report.FeatID = foundFeatID;
Report.FoundTime = foundTime;
Report.RTdelta = Report.RT-foundTime;
Report.FoundCompound = foundCompund;
Report.Score = foundScore;
