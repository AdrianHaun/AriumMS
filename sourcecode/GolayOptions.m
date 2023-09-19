classdef GolayOptions
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

    properties
        Container
        FrameSize
        Degree
        Flabel
        Dlabel
    end

    methods
        function obj = GolayOptions(CallingTab,CallingData)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Container = uipanel(CallingTab,...
                "Title","Golay Smoothing Parameters",...
                FontWeight="bold",...
                Position=[925,105,220,85],...
                TitlePosition="centertop",...
                BackgroundColor=[0.90,0.90,0.90],Enable="off");

            obj.FrameSize = uieditfield(obj.Container,"numeric",...
                Limits=[1,Inf],Value=CallingData.FrameSize,...
                Position=[125,30,90,20],...
                RoundFractionalValues="on",...
                ValueDisplayFormat='%.0f scans',...
                Tooltip="Moving Frame Size. Higher values smooth the signal more with an increase in computation time. Recommended setting: 15-30");
            obj.Flabel  = uilabel(obj.Container,"Text","Frame Size","Position",[5,30,110,20],"HorizontalAlignment","right");

            obj.Degree = uieditfield(obj.Container,"numeric",...
                Limits=[1,CallingData.FrameSize - 1],...
                Value=CallingData.Degree,...
                Position=[125,5,90,20],...
                RoundFractionalValues="on",...
                Tooltip="Specifies the degree of the polynomial fitted to the points in the moving frame. Recommended setting: 2-4");
            obj.Dlabel  = uilabel(obj.Container,"Text","Polynomial Degree","Position",[5,5,110,20],"HorizontalAlignment","right");
        end
    end

end
