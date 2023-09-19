classdef SampScaleOptions
    %UI objects for Sample and Group Scaling Options

    properties
        SelectButton
        ScalingLamp
        Container
        window
        InputTable
        Labels
        OKButton
        GroupValField
    end

    methods
        function obj = SampScaleOptions(CallingTab)
            %UNTITLED2 Construct an instance of this class
            obj.Container = uipanel(CallingTab,...
                "Title","Scaling Options",...
                FontWeight="bold",...
                Position=[10,10,220,65],...
                TitlePosition="centertop",...
                BackgroundColor=[0.90,0.90,0.90],Enable="off");

            obj.SelectButton = uibutton(obj.Container,"push",...
                "Position",[50,5,120,35],...
                Text="Define Factors");

            obj.ScalingLamp = uilamp(obj.Container,...
                "Color",[1,0,0],...
                Position=[200,47.5,15,15]);
        end
    end
end