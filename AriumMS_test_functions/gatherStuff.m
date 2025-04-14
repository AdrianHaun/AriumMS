rowName = ["ROI masses found","Features Found","Targets Found","FPR","average Mass Delta"];
ColNames = [app.Data.GroupName];

targetCells = cell(size(app.Data));
for n = 1:width(targetCells)
name = app.Data(n).GroupName;
vals = [targetDataTable.M,targetDataTable.(name)];
vals(vals(:,2)==0,:) = [];
targetCells{1,n} = vals;
end

stuff = zeros(5,width(app.Data));

for n = 1:width(app.Data)
    [responses,deltaMZ] = ObjectiveTest(app.Data(n).Output,targetCells{1,n},15);

    stuff(1,n) = numel(app.Data(n).ROIDataFileObj.ROImzVec);
    stuff(2,n) = responses(3);
    stuff(3,n) = responses(1);
    stuff(4,n) = responses(2);
    stuff(5,n) = mean(abs(deltaMZ),"omitmissing");
end

stuff = array2table(stuff,"RowNames",rowName,"VariableNames",ColNames)


clearvars -except targetDataTable app stuff