function [MSroi,time] = AlignPeaks(obj,MSroi,time,maxScan)
                % rearrange matrices
                [splitVar,~] = cellfun(@size,time);
                test= cellfun(@(x) sum(x~=0),time);
                [~,test] = max(test);
                TimeVec = time{test};
                MSroi=vertcat(MSroi{:});
                time=vertcat(time{:});
                [~,id]=max(MSroi);
                refTime=time(id);
                ShiftVal = [obj.maxshiftneg,obj.maxshiftpos];
                PW = obj.PulseWidth;
                WSR = obj.WindowSizeRatio;
                I = obj.Iterations;
                GS = obj.GridSteps;
                SS = obj.SearchSpace;
                parfor n=1:size(MSroi,2)
                    ROI=reshape(MSroi(:,n),maxScan,[]);
                    ROI=msalign(TimeVec,ROI,refTime(n),'MaxShift',ShiftVal,...
                        'WidthOfPulses',PW,'WindowSizeRatio',WSR,'Iterations',...
                        I,'GridSteps',GS,'SearchSpace',SS);
                    ROI(isnan(ROI))=0; %remove possible NaN
                    MStemp{1,n} = reshape(ROI,[],1);
                end
                MSroi = cell2mat(MStemp);
                MSroi = mat2cell(MSroi,splitVar);
                time = mat2cell(time,splitVar);
        end