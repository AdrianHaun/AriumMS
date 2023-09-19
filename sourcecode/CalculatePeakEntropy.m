function PeakEntropy = CalculatePeakEntropy(PeakData,EIC)
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
D=diff(EIC);
p=zeros(size(PeakData,1),1);
for n=1:size(PeakData,1)
    %extract peak range
    Peak = D(PeakData(n,2):PeakData(n,3));
    maxidx=PeakData(n,1)-PeakData(n,2);
    % check normal or variant point , variant point = 1
    premax=Peak(1:maxidx-1)<0;
    postmax=Peak(maxidx+1:end)>0;
    VarPoints=[premax; false; postmax];
    p(n,1)=sum(VarPoints)/numel(VarPoints);
end
PeakEntropy=-p.*log2(p)-(1-p).*log2(1-p);
PeakEntropy(isnan(PeakEntropy))=0;
end