function [MSroi_end,mzroi_end,time_end]=AutoROI(Peaklist,Timelist,mzerror,mzErrorUnit,minroi,thresh)
%%AutoROI Performs fully automated ROI search and augmentation.

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

%preallocate cell arrays
mzlist = cell(length(Peaklist),1);
MSroilist = cell(length(Peaklist),1);
%ROI search for every Sample
parfor d = 1 : length(Peaklist)
    P= Peaklist{d,1};
    T= Timelist{d,1};
    nrows=length(P);
    [mzlist{d,1},MSroilist{d,1},~]=ROIpeaks2(P,thresh,mzerror,mzErrorUnit,minroi,nrows,T);
end
if length(mzlist) == 1 %Skip Augmentation if only one Sample
    MSroi_end=MSroilist{1,1};
    mzroi_end=mzlist{1,1};
    time_end=Timelist{1,1};
else
    for i = 2:size(Peaklist,1)
        [MSroilist{1,1},mzlist{1,1},Timelist{1,1}] = MSroiaug2(MSroilist{1,1},MSroilist{i,1},mzlist{1,1},mzlist{i,1},mzerror,mzErrorUnit,thresh,Timelist{1,1},Timelist{i,1});
    end
    MSroi_end=MSroilist{1,1};
    mzroi_end=mzlist{1,1};
    time_end=Timelist{1,1};
end
end

