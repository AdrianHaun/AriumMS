function [MSroi,time] = AlignPeaks(obj,MSroi,time,maxScan)
% BSD 3-Clause License
% 
% Copyright (c) 2022, Adrian Haun
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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