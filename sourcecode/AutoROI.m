function [MSroi_end,mzroi_end,time_end]=AutoROI(Peaklist,Timelist,mzerror,minroi,thresh)
%%AutoROI Performs fully automated ROI search and augment.
%   Designd to be used in conjunction with ROIprocess
        %preallocate cell arrays
        mzlist = cell(length(Peaklist),1);
        MSroilist = cell(length(Peaklist),1);
        %ROI search for every Sample
            parfor d = 1 : length(Peaklist)
                P= Peaklist{d,1};
                T= Timelist{d,1};
                nrows=length(P);
                [mzroi,MSroi,~]=ROIpeaks(P,thresh,mzerror,minroi,nrows,T);
                mzlist{d,1} = mzroi;
                MSroilist{d,1} = MSroi;
            end
            if length(mzlist) == 1 %Skip Augmentation if only one Sample
                MSroi_end=MSroilist{1,1};
                mzroi_end=mzlist{1,1};
                time_end=Timelist{1,1};
            else
                ROIaugini
            end
            function ROIaugini %Augment first two sample ROIs
                MSroi1=MSroilist{1,1};
                MSroi2=MSroilist{2,1};
                mzroi1=mzlist{1,1};
                mzroi2=mzlist{2,1};
                time1=Timelist{1,1};
                time2=Timelist{2,1};
                [MSroi_aug,mzroi_aug,time_aug] = MSroiaug(MSroi1,MSroi2,mzroi1,mzroi2,mzerror,thresh,time1,time2);
                MS12=MSroi_aug;
                mz12=mzroi_aug;
                t12=time_aug;
                if length(Peaklist) == 2 %Skip further Augmentation if only 2 Samples
                    MSroi_end=MSroi_aug;
                    mzroi_end=mzroi_aug;
                    time_end=time_aug;
                else
                ROIAugloop
                end
                function ROIAugloop
                    for a = 1 : length(Peaklist)-2
                        i=a+2;
                        z=length(Peaklist);
                        mz=mzlist{i,1};
                        Ti=Timelist{i,1};
                        MS=MSroilist{i,1};
                        [MSroi_aug,mzroi_aug,time_aug] = MSroiaug(MS12,MS,mz12,mz,mzerror,thresh,t12,Ti);                 
                        if z>i %Check if loop must continue
                            MS12=MSroi_aug; 
                            mz12=mzroi_aug;
                            t12=time_aug;
                        else
                            MSroi_end=MSroi_aug;
                            mzroi_end=mzroi_aug;
                            time_end=time_aug;
                        end
                    end
                end
            end
end

