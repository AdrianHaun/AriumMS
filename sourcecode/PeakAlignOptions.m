classdef PeakAlignOptions
%% Object for storing Peak Alignment UI elements
    properties
        Container 
        MaxShiftFieldPos
        MaxShiftFieldNeg
        PulseWidthField
        WindowSizeField
        SearchSpaceDrop
        IterationField
        GridStepField
        Labels
    end

    methods
        function obj = PeakAlignOptions(CallingTab,CallingData)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here            
            obj.Container = uipanel(CallingTab,"Title","Peak Alignment Parameters", ...
                "BackgroundColor",[0.90,0.90,0.90], ...
                "Position",[240,10,220,200], ...
                "TitlePosition","centertop",...
                "FontWeight","bold",Enable="off");

            obj.MaxShiftFieldNeg = uieditfield(obj.Container,"numeric", ...
                "Position",[125,155,90,20], ...
                "Limits",[-Inf,0], ...
                "ValueDisplayFormat","%.0f scans",...
                "Value",CallingData.maxshiftneg, ...
                "RoundFractionalValues","on",...
                Tooltip="Maximum allowed shift of peaks towards lower retention times.");
            obj.MaxShiftFieldPos = uieditfield(obj.Container,"numeric", ...
                "Position",[125,130,90,20], ...
                "Limits",[1,Inf], ...
                "ValueDisplayFormat","%.0f scans",...
                "Value",CallingData.maxshiftpos, ...
                "RoundFractionalValues","on",...
                Tooltip="Maximum allowed shift of peaks towards higher retention times.");
            obj.Labels(1) = uilabel(obj.Container,"Text","Maximum Shift","Position",[5,142.5,110,20],HorizontalAlignment="right");

            obj.PulseWidthField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,105,90,20], ...
                "Limits",[1,Inf], ...
                "ValueDisplayFormat","%.0f scans",...
                "Value",CallingData.PulseWidth, ...
                "RoundFractionalValues","on",...
                Tooltip="Width of the synthetic gauss peak at 60.65% of it's maximum height");
            obj.Labels(2) = uilabel(obj.Container,"Text","Pulse Width","Position",[5,105,110,20],HorizontalAlignment="right");

            obj.WindowSizeField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,80,90,20], ...
                "Limits",[0.01,Inf], ...
                "ValueDisplayFormat","%11.4g",...
                "Value",CallingData.WindowSizeRatio, ...
                "RoundFractionalValues","off",...
                Tooltip="Scaling factor that determines the size of the window around every alignment peak. The synthetic signal is compared to the input signal only within these regions, which saves computation time.");
            obj.Labels(3) = uilabel(obj.Container,"Text","Window Size Ratio","Position",[5,80,110,20],HorizontalAlignment="right");
            
            obj.SearchSpaceDrop = uidropdown(obj.Container,...
                "Value",CallingData.SearchSpace,...
                "Items",["Evenly Spaced Lattice","Random Latin Hypercube"],...
                "ItemsData",["regular","latin"],...
                "Position",[125,55,90,20],...
                Tooltip="Search space type");
            obj.Labels(4) = uilabel(obj.Container,"Text","Search Space","Position",[5,55,110,20],HorizontalAlignment="right");

            obj.IterationField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,30,90,20], ...
                "Limits",[1,Inf], ...
                "ValueDisplayFormat","%.0f",...
                "Value",CallingData.Iterations, ...
                "RoundFractionalValues","on",...
                Tooltip="Number of refining iterations. At every iteration, the search grid is scaled down to improve the estimates.");
            obj.Labels(5) = uilabel(obj.Container,"Text","Iterations","Position",[5,30,110,20],HorizontalAlignment="right");

            obj.GridStepField = uieditfield(obj.Container,"numeric", ...
                "Position",[125,5,90,20], ...
                "Limits",[1,Inf], ...
                "ValueDisplayFormat","%.0f",...
                "Value",CallingData.GridSteps, ...
                "RoundFractionalValues","on",...
                Tooltip="Specifies the number of steps for the search grid. At every iteration, the search area is divided by GridSteps^2.");
            obj.Labels(6) = uilabel(obj.Container,"Text","Grid Steps","Position",[5,5,110,20],HorizontalAlignment="right");
        end
    end
end