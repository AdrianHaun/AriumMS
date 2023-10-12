classdef GolayOptions
%% Object for storing Golay Smoothing Option UI elements

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
