function valuesFiltered = AdductFilterAlgo(values,AdductRules,mztol,CosSim)
%% AdductFilterAlgo Filters Adduct Peaks from Internal AriumMS integration results
%
%   Calculates possible non Adduct (Base) m/z for each extracted m/z,
%   then finds matching masses in original list. Peak lists of Adduct and
%   Base m/z are then compared to have peaks with the same RT (within
%   2 sec). Matching peaks shapes are compared using Cosine Similarity
%   (>=0.85 default)

%Preparation
RTs=cellfun(@(X) X(:,2),values(4,:),'UniformOutput',false); %extract Retentiontimes
ranges=cellfun(@(X) X(:,2:3),values(3,:),'UniformOutput',false); %extract Peak ranges
mz=cell2mat(values(1,:));

%build possible Adduct mass list
Case=zeros(length(AdductRules),length(mz));
for n=1:length(AdductRules)
    Case(n,:)=AdductRules{n}(mz);
end

%compare case masses to initial mz list
[~,Case]=ismembertol(Case,mz,mztol,'DataScale',1,'OutputAllIndices',true);

%%Check Adducts
%transform cell array to array
rows=max(sum(cellfun(@numel,Case)));
idxMat=zeros(rows,length(mz));
for k=1:length(mz)
    val=vertcat(Case{:,k});
    idxMat(1:length(val),k)=val;
end
%find matching RT indices
matchRT=cell(size(idxMat));
for n=1:size(idxMat,2)
    for k=1:size(idxMat,1)
        if idxMat(k,n)~=0
            %check RT
            [~,matchRT{k,n}]=ismembertol(RTs{n},RTs{idxMat(k,n)},2,'DataScale',1);
        end
    end
end
%remove lists with only zeros
check=cellfun(@(X) all(X(:)==0), matchRT);
matchRT(check)={[]};
idxMat(check)=0;

%Preallocate storage variable of identified Adduct variables
found=cell(1,length(mz));

for n=1:size(matchRT,2) %Adduct mz loop
    test=vertcat(matchRT{:,n}); 
    if isempty(test) == false       %check if any matching Peaks present, if true extract relevant Adduct data
        AdductEIC = values{6,n};    
        AdductRanges=ranges{1,n};   
        for k=1:size(matchRT,1)     % possible Base mz loop
            TestSet=matchRT{k,n};
            if isempty(TestSet) == false            % check if matching Peak present, if true extract relevant data
                BaseEIC=values{6,idxMat(k,n)};      %base m/z EIC
                BaseRanges=ranges{1,idxMat(k,n)};   %base m/z peak ranges
                CheckList=matchRT{k,n};             %List of Peaks that match Adduct rt
                %preallocate variables
                idx=false(size(CheckList,1),1);
                isAddukt=false(size(idx,1),size(matchRT,1));
                for p=1:size(CheckList,1)       % loop over matching Peaks
                    if CheckList(p,1)~=0        % check if Peak in Basepeak List has matching Peak in Adduct List
                        RangeA=AdductRanges(p,:);
                        RangeB=BaseRanges(CheckList(p,1),:);
                        EICA=AdductEIC;
                        EICB=BaseEIC;
                        %Set the intensities outside the considered range to 0
                        EICA(1:RangeA(1)-1)=0;
                        EICA(RangeA(2)+1:end)=0;
                        EICB(1:RangeB(1)-1)=0;
                        EICB(RangeB(2)+1:end)=0;
                        idx(p)=sum(EICB.*EICA)/(sqrt(sum(EICB.^2))*sqrt(sum(EICA.^2)))>=CosSim;
                    end
                end
                isAddukt(:,k)=idx;
            end
        end
        found{1,n}=any(isAddukt,2);
    end
end
%remove identified Adducts
valuesFiltered=Removify(values,found);

%%support function
    function values=Removify(values,found) %removes confirmed Adduct from Input list and stores number of removed peaks
        for i=1:length(found)
            if isempty(found{i})==false
                values{13,i}=sum(found{i});
                values{2,i}(found{1,i},:)=[];
                values{3,i}(found{1,i},:)=[];
                values{4,i}(found{1,i},:)=[];
                values{5,i}(found{1,i},:)=[];
            else
                values{13,i}=0;
            end
        end
    end
end
