classdef BaseCorrOptions
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
        WindowSizeField
        StepSizeField
        RegressionDrop
        EstimationDrop
        SmoothDrop
        QuantilValField
        Labels
    end

    methods
        function obj = BaseCorrOptions(CallingTab,CallingData)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here            
            obj.Container = uipanel(CallingTab,"Title","Baseline Correction Parameters", ...
                "BackgroundColor",[0.90,0.90,0.90], ...
                "Position",[465,10,220,180], ...
                "TitlePosition","centertop",...
                FontWeight="bold",Enable="off");

            obj.WindowSizeField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,130,90,20], ...
                "Limits",[1,Inf], ...
                "ValueDisplayFormat","%.0f scans",...
                "Value",CallingData.WindowSize, ...
                "RoundFractionalValues","on",...
                Tooltip="Window size for baseline estimation. Recommended setting: 150-250");
            obj.Labels(1) = uilabel(obj.Container,"Text","Window Size","Position",[5,130,110,20],"HorizontalAlignment","right");

            obj.StepSizeField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,105,90,20], ...
                "Limits",[1,Inf], ...
                "ValueDisplayFormat","%.0f scans",...
                "Value",CallingData.StepSize, ...
                "RoundFractionalValues","on",...
                Tooltip="Defines placement of baseline estimation windows every X scans. Recommended setting: Shifting window size or less");
            obj.Labels(2) = uilabel(obj.Container,"Text","Step Size","Position",[5,105,110,20],"HorizontalAlignment","right");
            
            obj.RegressionDrop = uidropdown(obj.Container,...
                "Value",CallingData.RegressionMethod,...
                "Items",["Shape-preserving piecewise cubic interpolation","Linear Interpolation","Spline Interpolation"],...
                "ItemsData",["pchip","linear","spline"],...
                "Position",[125,80,90,20],...
                Tooltip="Method to regress the window estimated points to a soft curve, specified as one of the following: Shape-preserving piecewise cubic interpolation, linear interpolation, Spline interpolation. Recommended setting: Shape-preserving piecewise cubic interpolation");
            obj.Labels(3) = uilabel(obj.Container,"Text","Regression Method","Position",[5,80,110,20],"HorizontalAlignment","right");

            obj.EstimationDrop = uidropdown(obj.Container,...
                "Value",CallingData.EstimationMethod,...
                "Items",["Quantile","Expectation-Maximization"],...
                "ItemsData",["quantile","em"],...
                "Position",[125,55,90,20],...
                Tooltip="Method to find likely baseline value in every window, specified as one of the following: Quantile, Expectation-Maximization: Every sample is the independent and identically distributed draw of any of two normal distributed classes (background or peaks). Because the class label is hidden, the distributions are estimated with an Expectation-Maximization algorithm. The ultimate baseline value is the mean of the background class");
            obj.Labels(4) = uilabel(obj.Container,"Text","Estimation Method","Position",[5,55,110,20],"HorizontalAlignment","right");

            obj.SmoothDrop = uidropdown(obj.Container,...
                "Value",CallingData.SmoothMethod,...
                "Items",["No Smoothing","Linear fit","Quadratic fit","Robust linear fit","Robust quadratic fit"],...
                "ItemsData",["none","lowess","loess","rlowess","rloess"],...
                "Position",[125,30,90,20],...
                Tooltip="Method to smooth the curve of estimated points. Recomemded setting: none");
            obj.Labels(5) = uilabel(obj.Container,"Text","Smoothing Method","Position",[5,30,110,20],"HorizontalAlignment","right");

            obj.QuantilValField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,5,90,20], ...
                "Limits",[0,1], ...
                "ValueDisplayFormat","%11.4g",...
                "Value",CallingData.QuantilVal, ...
                "RoundFractionalValues","off",...
                Tooltip="Quantile value for Quantile estimation method. Recommended setting: 0.1");
            obj.Labels(6) = uilabel(obj.Container,"Text","Quantile Value","Position",[5,5,110,20],"HorizontalAlignment","right");

        end
    end
end