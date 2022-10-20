function valuesFiltered = IsotopeFilterAlgo(values,IsotopeRules,mztol)
%% IsotopeFilterAlgo Filters Isotope Peaks from Internal AriumMS integration results
%
%   Calculates possible non Isotope (Base) m/z forch each extracted m/z,
%   then finds matching masses in original list. Peak lists of Isotope and
%   Base m/z are then compared to have peaks with the same RT (within
%   2 sec). Matching peaks shapes are compared using Cosine Similarity
%   (>=0.85), and Base m/z intensity adjusted by the relative Isotope
%   occurence must be within 10% of the Isotope intensity

%Preparation
RTs=cellfun(@(X) X(:,2),values(4,:),'UniformOutput',false); %extract Retentiontimes
ranges=cellfun(@(X) X(:,2:3),values(3,:),'UniformOutput',false); %extract Peak ranges
mz=cell2mat(values(1,:));
Case=mz-IsotopeRules(:,1); %build possible Isotope mass list
%compare case masses to initial mz list
[~,Case]=ismembertol(Case,mz,mztol,'DataScale',1,'OutputAllIndices',true);

%build natural isotope occurence factor list 
IsoOCC=cell(size(Case));
for n=1:size(Case,1)
    IsoOCC(n,:)=cellfun(@(X) ones(numel(X),1)*IsotopeRules(n,2), Case(n,:), 'UniformOutput', false);
end
%%Check Adducts
%transform cell array to array
rows=max(sum(cellfun(@numel,Case)));
idxMat=zeros(rows,length(mz),2);
for k=1:length(mz)
    val=vertcat(Case{:,k});
    idxMat(1:length(val),k,1)=val;
    val=vertcat(IsoOCC{:,k});
    idxMat(1:length(val),k,2)=val;
end
%find matching RT indices
matchRT=cell(size(idxMat(:,:,1)));
for n=1:size(idxMat,2)
    for k=1:size(idxMat,1)
        if idxMat(k,n)~=0
            %check RT
            [~,matchRT{k,n}]=ismembertol(RTs{n},RTs{idxMat(k,n,1)},2,'DataScale',1);
        end
    end
end
%remove lists with only zeros
check=cellfun(@(X) all(X(:)==0), matchRT);
matchRT(check)={[]};
check=repmat(check,1,1,2);
idxMat(check)=0;

%Preallocate storage variable of identified Isotopes variables
found=cell(1,length(mz));

for n=1:size(matchRT,2) %Isotope mz loop
    test=vertcat(matchRT{:,n}); 
    if isempty(test) == false %check if any matching Peaks present, if true extract relevant Isotope data
        IsotopeEIC = values{6,n};
        IsotopeRanges=ranges{1,n};
        for k=1:size(matchRT,1) % possible Base mz loop
            TestSet=matchRT{k,n};
            if isempty(TestSet) == false % check if matching Peak present, if true extract relevant data
                BaseEIC=values{6,idxMat(k,n,1)};    %base m/z EIC
                BaseRanges=ranges{1,idxMat(k,n,1)}; %Base m/z peak ranges
                CheckList=matchRT{k,n};             %List of Peaks that match Isotope rt
                factor=idxMat(k,n,2);               %relative Intensity of Isotopes
                %preallocate variables
                idx=false(size(CheckList,1),1);
                isISO=false(size(idx,1),size(matchRT,1));
                for p=1:size(CheckList,1)   % loop over matching Peaks
                    if CheckList(p,1)~=0 % check if Peak in Basepeak List has matching Peak in Isotope List
                        RangeA=IsotopeRanges(p,:);
                        RangeB=BaseRanges(CheckList(p,1),:);
                        EICA=IsotopeEIC;
                        EICB=BaseEIC;
                        %Set the intensities outside the considered range to 0
                        EICA(1:RangeA(1)-1)=0;
                        EICA(RangeA(2)+1:end)=0;
                        EICB(1:RangeB(1)-1)=0;
                        EICB(RangeB(2)+1:end)=0;
                        tolA= max(EICA)*+max(EICA)*0.10;
                        tolB= max(EICA)*-max(EICA)*0.10;
                        idx(p)=sum(EICB.*EICA)/(sqrt(sum(EICB.^2))*sqrt(sum(EICA.^2)))>=0.85 & max(EICB)*factor<=tolA & max(EICB)*factor>=tolB; %desision Cosine Similarity and  mainPeak Intensity adjusted with relative Isotope Occurence is inside 20% of isotopeintensity
                    end
                end
                isISO(:,k)=idx;
            end
        end
        found{1,n}=any(isISO,2);
    end
end
%remove identified Isotopes
valuesFiltered=Removify(values,found);

%%support function

    function values=Removify(values,found) %removes confirmed Isotope from Input list and stores number of removed peaks
        for i=1:length(found)
            if isempty(found{i})==false
                values{12,i}=sum(found{i});
                values{2,i}(found{1,i},:)=[];
                values{3,i}(found{1,i},:)=[];
                values{4,i}(found{1,i},:)=[];
                values{5,i}(found{1,i},:)=[];
            else
                values{12,i}=0;
            end
        end
    end
end
