function [MSroi,time] = SmoothPeaks(obj,MSroi,time)
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
Frame = obj.FrameSize;
Deg = obj.Degree;
parfor id = 1:size(MSroi,1)
    %depad Array
    oldSize=size(MSroi{id});
    MSroiTemp = MSroi{id};
    [MSroiTemp,timeTemp] = depadArrays(MSroiTemp,time{id});
    MSroiTemp = mssgolay(timeTemp,MSroiTemp,'Span',Frame,'Degree',Deg);
    %remove negativ, NaN and repad Array
    MSroiTemp=max(MSroiTemp,0);
    MSroiTemp(isnan(MSroiTemp))=0;
    [MSroi{id},time{id}] = repadArrays(MSroiTemp,timeTemp,oldSize);
    % set possible negative values to 0
    MSroi{id} = max(MSroi{id},0);
end

end
function [ArrayOut,VecOut] = depadArrays(ArrayIn,VecIn)
idx = VecIn~=0;
ArrayOut = ArrayIn(idx,:);
VecOut = VecIn(idx);
end

function [ArrayOut,VecOut] = repadArrays(ArrayIn,VecIn,PadSize)
OrgSize=size(ArrayIn);
ArrayOut = padarray(ArrayIn,PadSize-OrgSize,0,'post');
VecOut = padarray(VecIn,PadSize(1)-OrgSize(1),0,'post');
end