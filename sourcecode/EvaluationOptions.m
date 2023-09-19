classdef EvaluationOptions
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
        Container                       matlab.ui.container.Panel
        DataVisualizationSubPanel       matlab.ui.container.Panel
        MultivariateAnalysisSubPanel    matlab.ui.container.Panel
        TICButton                       matlab.ui.control.Button
        ROIButton                       matlab.ui.control.Button
        IntegrationButton               matlab.ui.control.Button
        EntropyButton                   matlab.ui.control.Button
        SNButton                        matlab.ui.control.Button
        RejectedButton                  matlab.ui.control.Button
        ScatterButton                   matlab.ui.control.Button
        DistributionButton              matlab.ui.control.Button
        pDistributionButton             matlab.ui.control.Button
        CloudButton                     matlab.ui.control.Button
        ProbabilityButton               matlab.ui.control.Button
        VolcanoButton                   matlab.ui.control.Button
        HeatmapButton                   matlab.ui.control.Button
        PCAButton                       matlab.ui.control.Button
        ExportButton                    matlab.ui.control.Button
        PlotSpace                       matlab.graphics.layout.TiledChartLayout
        PlotTiles                       matlab.graphics.axis.Axes
    end
    
    methods
        function obj = EvaluationOptions(CallingApp)
            %EVALUATIONOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Container = uipanel(CallingApp.MainWindow,"Title","EvaluationOptions","BackgroundColor",[0.94,0.94,0.94],"Position",[10,10,1155,655]);
            obj.DataVisualizationSubPanel = uipanel(obj.Container,"Title","Data Visualization","BackgroundColor",[0.9,0.9,0.9],"Position",[10,10,190,615]);
            obj.MultivariateAnalysisSubPanel = uipanel(obj.Container,"Title","Data Analysis","BackgroundColor",[0.9,0.9,0.9],"Position",[210,10,550,85]);
           %% Visualization
            obj.TICButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","Plot TIC",...
                "Position",[10,540,170,45],Enable="off");
            obj.ROIButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","Plot found ROIs",...
                "Position",[10,490,170,45],Enable="off");
            obj.IntegrationButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","Plot Integration",...
                "Position",[10,440,170,45],Enable="off");
            obj.EntropyButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text",["Peak Entropy", "Distribution"],...
                "Position",[10,390,170,45],Enable="off");
            obj.SNButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text",["Signal to Noise"; "Distribution"],...
                "Position",[10,340,170,45],Enable="off");
            obj.RejectedButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text",["Rejected Peaks";"by Filter"],...
                "Position",[10,290,170,45],Enable="off");
            obj.ScatterButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","Scatter Plot",...
                "Position",[10,240,170,45],Enable="off");
            obj.CloudButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","Cloud Plot",...
                "Position",[10,190,170,45],Enable="off");
            obj.ProbabilityButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","Probability Plot",...
                "Position",[10,140,170,45],Enable="off");
             obj.DistributionButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","Feature Distribution",...
                "Position",[10,90,170,45],Enable="off");
              obj.pDistributionButton = uibutton(obj.DataVisualizationSubPanel,...
                "Text","p-Value Distribution",...
                "Position",[10,40,170,45],Enable="off");
            %% Data Analysis
            obj.VolcanoButton = uibutton(obj.MultivariateAnalysisSubPanel,...
                "Text","Volcano Plot",...
                "Position",[10,10,170,45],Enable="off");
             obj.HeatmapButton = uibutton(obj.MultivariateAnalysisSubPanel,...
                "Text","Heatmap",...
                "Position",[190,10,170,45],Enable="off");
             obj.PCAButton = uibutton(obj.MultivariateAnalysisSubPanel,...
                "Text","PCA",...
                "Position",[370,10,170,45],Enable="off");
            %% Export and Info
            obj.ExportButton = uibutton(obj.Container,...
                "Text","Export Results",...
                "Position",[975,10,170,85],FontSize=15,FontWeight="bold",Enable="off");
            
            %% Plot Area
            numRows = ceil(size(CallingApp.Data,2)/2);
            obj.PlotSpace = tiledlayout(obj.Container,2,numRows,"TileSpacing","tight",Units="pixels",Position=[255 145 880 450],Visible="on");
        end
    end
end

