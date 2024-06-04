function SortedScores = SortBetweenScores(GroupScores)

Groups = unique(vertcat(GroupScores{1,:}));
SortedScores = cell(1,size(Groups,1));

for n = 1:size(Groups,1)
    GroupID = Groups(n);
    id = cellfun(@(X) any(X==GroupID,2),GroupScores(1,:));
    Scores = GroupScores(:,id);
    TempScores = [];
    for i = size(Scores,2)
        LocalScores = Scores{2,i};
        LocalID = Scores{1,i} == GroupID;
        %get column id
        LocalID = [true,LocalID];
        TempScores = [TempScores;LocalScores(:,LocalID)];
    end
    SortedScores{1,n} = TempScores;
end