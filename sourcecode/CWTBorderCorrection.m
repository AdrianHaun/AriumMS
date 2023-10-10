function CorrectedPeakData = CWTBorderCorrection(PeakData,EIC,smoothedEIC)
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

EIC=full(EIC);
smoothedEIC = full(smoothedEIC);
noise=EIC-smoothedEIC;
%friction border correction
smoothedEIC=rescale(smoothedEIC,0,1);
IntThresh=0.0001;
%normalize thresh
IntThresh=(max(smoothedEIC)-min(smoothedEIC))*IntThresh;
for n=1:size(PeakData,1)
    u=PeakData(n,3);
    while u<numel(smoothedEIC)-1 %upper bond
        if smoothedEIC(u)-smoothedEIC(u+1)>IntThresh
            u=u+1;
        else
            PeakData(n,3)=u;
            break
        end
    end
    l=PeakData(n,2);
    while l>1 %lower bond
        if smoothedEIC(l)-smoothedEIC(l-1)>IntThresh
            l=l-1;
        else
            PeakData(n,2)=l;
            break
        end
    end
end
%moving STD correction to remove errors from smoothed peaks
for n=1:size(PeakData,1)
    u=PeakData(n,3);
    while u>1 && u<numel(EIC)-8 %upper bond
        if abs(EIC(u-1)-mean(EIC(u-1:u+7,1))) <= std(noise(u-1:u+7,1)) || round(EIC(u-1,1)-EIC(u,1))==0
            u=u-1;
        else
            PeakData(n,3)=u;
            break
        end
    end
    l=PeakData(n,2);
    while l<=numel(EIC)-1 && l>8 %lower bond
        if abs(EIC(l+1,1)-mean(EIC(l-7:l+1,1))) <= std(noise(l-7:l+1,1)) || round(EIC(l,1)-EIC(l+1,1))==0
            l=l+1;
        else
            PeakData(n,2)=l;
            break
        end
    end
end
CorrectedPeakData = PeakData;
end