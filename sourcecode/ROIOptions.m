classdef ROIOptions
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Container       
        ThreshField     matlab.ui.control.NumericEditField
        mzerrorField    matlab.ui.control.NumericEditField
        minROIField     matlab.ui.control.NumericEditField
        StartTimeField  matlab.ui.control.NumericEditField
        EndTimeField    matlab.ui.control.NumericEditField
        Labels
    end

    methods
        function obj = ROIOptions(CallingTab,CallingData)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
         
            obj.Container = uipanel(CallingTab,"Title","ROI Parameters", ...
                "BackgroundColor",[0.90,0.90,0.90], ...
                "Position",[10,80,220,150], ...
                "FontWeight","bold",...
                "TitlePosition","centertop",Enable="off");
            
            obj.ThreshField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,105,90,20], ...
                "Limits",[1,Inf], ...
                "ValueDisplayFormat","%.0f",...
                "Value",CallingData.thresh, ...
                "RoundFractionalValues","on");
            obj.Labels(1) = uilabel(obj.Container,"Text","Intensity Threshold","Position",[5,105,110,20],HorizontalAlignment="right");

            obj.mzerrorField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,80,90,20], ...
                "Limits",[0,Inf], ...
                "LowerLimitInclusive","off",...
                "ValueDisplayFormat","%.4f Da",...
                "Value",CallingData.mzerror);
            obj.Labels(2) = uilabel(obj.Container,"Text","MZ error","Position",[5,80,110,20],HorizontalAlignment="right");

            obj.minROIField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,55,90,20], ...
                "Limits",[1,Inf], ...
                "LowerLimitInclusive","off",...
                "ValueDisplayFormat","%.0f scans",...
                "Value",CallingData.minroi, ...
                "RoundFractionalValues","on");
            obj.Labels(3) = uilabel(obj.Container,"Text","Minimum ROI size","Position",[5,55,110,20],HorizontalAlignment="right");

            obj.StartTimeField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,30,90,20], ...
                "Limits",[CallingData.Start, CallingData.End], ...
                "LowerLimitInclusive","on",...
                "ValueDisplayFormat","%.1f s",...
                "Value",CallingData.Start);
            obj.Labels(4) = uilabel(obj.Container,"Text","Start Time","Position",[5,30,110,20],HorizontalAlignment="right");

            obj.EndTimeField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,5,90,20], ...
                "Limits",[CallingData.Start + 0.1 ,CallingData.End], ...
                "UpperLimitInclusive","on",...
                "ValueDisplayFormat","%.1f s",...
                "Value",CallingData.End);
            obj.Labels(5) = uilabel(obj.Container,"Text","End Time","Position",[5,5,110,20],HorizontalAlignment="right");
        end        
    end
end