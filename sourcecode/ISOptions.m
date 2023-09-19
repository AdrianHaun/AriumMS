classdef ISOptions
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
        NumberSelect
        MassCalCheck
        ApplyDrop
        ApplyLabel matlab.ui.control.Label
        OrderDrop
        OrderLabel matlab.ui.control.Label
        DefineButton
        StatusLamp
        Labels
        window
        InputTable
        OKButton
        oldDat
    end

    methods
        function obj = ISOptions(CallingTab,CallingData)
            obj.oldDat = CallingData.ISDat;
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Container = uipanel(CallingTab,...
                "Title","Internal Standard Options",...
                FontWeight="bold",...
                Position=[695,10,220,180],...
                TitlePosition="centertop",...
                BackgroundColor=[0.90,0.90,0.90],Enable="off");

            obj.NumberSelect = uispinner(obj.Container,...
                Limits=[1,Inf],Value=CallingData.numISTD,...
                Step=1,...
                Position=[125,125,90,20],...
                ValueDisplayFormat="%.0f",...
                RoundFractionalValues="on",...
                Tooltip="Define the number of internal standards. If multiple standards are selected, feature intensities are normalized with the internal standard with the closest retention time.");
            obj.Labels(1)  = uilabel(obj.Container,"Text",["Number of","Internal Standards"],"Position",[5,115,110,40],"HorizontalAlignment","right",WordWrap=true);
            
            obj.MassCalCheck = uicheckbox(obj.Container,...
                Value=CallingData.MassCal,...
                Position=[10,25,175,20],...
                Text="Mass Correction",...
                Tooltip="Corrects the feature masses with the difference of the specified IS mass and the found IS mass. Feature masses are corrected with the internal standard with the closest mass.");

            obj.ApplyDrop = uidropdown(obj.Container,...
                "Value",CallingData.ISApply,...
                "Items",["Samples and Blanks","Samples only"],...
                "ItemsData",["S&B","SOnly"],...
                "Position",[105,95,110,20],"Enable",CallingData.BLKSubtraction,...
                Tooltip="Specify which files to apply IS normalization to. Use ""Samples only"" if the blank files do not contain IS.");
            obj.ApplyLabel = uilabel(obj.Container,"Text","Apply to","Position",[5,95,90,20],"Enable",CallingData.BLKSubtraction,HorizontalAlignment="right");

            obj.OrderDrop = uidropdown(obj.Container,...
                "Value",CallingData.ISOrder,...
                "Items",["Blank then IS","IS then blank"],...
                "ItemsData",["BlankIS","ISBlank"],...
                "Position",[105,65,110,20],...
                "Enable",CallingData.BLKSubtraction,...
                Tooltip="Specify the order of internal standard correction and blank subtraction.");
            obj.OrderLabel = uilabel(obj.Container,"Text","Normalization order","Position",[5,55,90,40],"Enable",CallingData.BLKSubtraction,HorizontalAlignment="right",WordWrap="on");

            obj.DefineButton = uibutton(obj.Container,...
                "Text","Define ISTD",...
                "Position",[125,15,90,40]);

            obj.StatusLamp = uilamp(obj.Container,...
                "Position",[200,162.5,15,15]);
            if isempty(obj.oldDat) == true
               obj.StatusLamp.Color=[1.00,0.00,0.00];
            else
                obj.StatusLamp.Color=[0.39,0.83,0.07];
            end
        end
        
    end
end