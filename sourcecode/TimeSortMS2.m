function SortedSpectra = TimeSortMS2(SpectraCells,Times)
    [~,index] = sort(Times,'ascend');
    SortedSpectra = SpectraCells(index);
end