classdef ISOptions
    %% Object for storing Internal Standard option UI elements

    properties
        Container
        NumberSelect
        MassCalCheck
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

            if nargin > 0
                
                obj.oldDat = CallingData.internalStandardData;

                obj.Container = uipanel(CallingTab,...
                    "Title","Internal Standard Options",...
                    FontWeight="bold",...
                    Position=[695,10,220,180],...
                    TitlePosition="centertop",...
                    BackgroundColor=[0.90,0.90,0.90],Enable="off");

                obj.NumberSelect = uispinner(obj.Container,...
                    Limits=[1,Inf],Value=CallingData.nInternalStandard,...
                    Step=1,...
                    Position=[125,125,90,20],...
                    ValueDisplayFormat="%.0f",...
                    RoundFractionalValues="on",...
                    Tooltip="Define the number of internal standards. If multiple standards are selected, feature intensities are normalized with the internal standard with the closest retention time.");
                obj.Labels(1)  = uilabel(obj.Container,"Text",["Number of","Internal Standards"],"Position",[5,115,110,40],"HorizontalAlignment","right",WordWrap=true);

                obj.MassCalCheck = uicheckbox(obj.Container,...
                    Value=CallingData.useISMassCorrection,...
                    Position=[10,25,175,20],...
                    Text="Mass Correction",...
                    Tooltip="Corrects the feature masses with the difference of the specified IS mass and the found IS mass. Feature masses are corrected with the internal standard with the closest mass.");

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
end