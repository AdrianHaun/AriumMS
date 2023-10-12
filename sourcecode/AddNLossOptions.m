classdef AddNLossOptions
%% Object for storing Adduct and Neutral loss UI elements
    properties
        SelectButton
        AddNLossLamp
        Container
        window
    end

    methods
        function obj = AddNLossOptions(CallingTab)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Container = uipanel(CallingTab,...
                "Title","   Adduct & Neutral Loss Options",...
                FontWeight="bold",...
                Position=[925,10,220,80],...
                TitlePosition="lefttop",...
                BackgroundColor=[0.90,0.90,0.90],Enable="off");

            obj.SelectButton = uibutton(obj.Container,"push",...
                "Position",[50,10,120,40],...
                Text="Define Filter");

            obj.AddNLossLamp = uilamp(obj.Container,...
                "Color",[1,0,0],...
                Position=[200,62.5,15,15]);
        end
    end
end