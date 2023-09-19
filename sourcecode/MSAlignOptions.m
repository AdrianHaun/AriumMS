classdef MSAlignOptions
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
        QuantilField
        EstimationDrop
        CorrectionDrop
        Labels
    end

    methods
        function obj = MSAlignOptions(CallingTab,CallingData)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.Container = uipanel(CallingTab,"Title","Spectral Alignment Parameters", ...
                "BackgroundColor",[0.90,0.90,0.90], ...
                "Position",[240,215,220,105], ...
                "TitlePosition","centertop",...
                "FontWeight","bold",Enable="off");

            obj.QuantilField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,55,90,20], ...
                "Limits",[0,1], ...
                "ValueDisplayFormat","%11.4g",...
                "Value",CallingData.mzQuantil, ...
                "RoundFractionalValues","off",...
                Tooltip="Determines which peaks are selected by the estimation method to create the vector of common m/z values (CMZ).");
            
            obj.Labels(1) = uilabel(obj.Container,"Text","Quantil Value","Position",[5,55,110,20],"HorizontalAlignment","right");

            obj.EstimationDrop = uidropdown(obj.Container, ...
                "Items",["Histogram","Regression"],...
                ItemsData=["histogram","regression"], ...
                Value=CallingData.mzEstimMethod,...
                Position=[125,30,90,20],...
                Tooltip="Specifies the method used to estimate the CMZ values. Choises are: Histogram: Peak locations are clustered using a kernel density estimation approach. The peak ion intensity is used as a weighting factor. The center of all the clusters conform to the CMZ vector; Regression  Takes a sample of the distances between observed significant peaks and regresses the inter-peak distance to create the CMZ vector with similar inter-element distances.");
            
            obj.Labels(2) = uilabel(obj.Container,"Text","Estimation Method","Position",[5,30,110,20],"HorizontalAlignment","right");

            obj.CorrectionDrop = uidropdown(obj.Container, ...
                "Items",["Nearest Neighbor","Shortest Path"],...
                ItemsData=["nearest-neighbor","shortest-path"], ...
                Value=CallingData.mzCorrectionMethod,...
                Position=[125,5,90,20],...
                Tooltip="Specifies the method to align each mass spectra. Choises are: Nearest-Neighbor: For each common peak in the CMZ vector, its counterpart in each peak list is the peak that is closest to the common peak's m/z value; Shortest-Path: For each common peak in the CMZ vector, its counterpart in each peak list is selected using the shortest path algorithm.");
            
            obj.Labels(3) = uilabel(obj.Container,"Text","Correction Method","Position",[5,5,110,20],"HorizontalAlignment","right");
        end
    end
end