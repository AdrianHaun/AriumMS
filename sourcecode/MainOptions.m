classdef MainOptions
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Container
        PolaritySwitch
        PolLabel
        BlankTickBox
        PeakAlignmentTickBox
        MSAlignmentTickBox
        BaseCorrTickBox
        GolayTickBox
        ScalingTickBox
        ISTDTickBox
        IsotopeFilterTickBox
        AddNLossTickBox
        ContaminantTickBox
        Labels
    end

    methods
        function  obj = MainOptions(CallingTab,CallingData)
            % Construct an instance of this class
            %   Construct Main Options Panel for new AriumMS Tab
            obj.Container = uipanel(CallingTab,"Title","Main Options",...
                FontWeight="bold",...
                TitlePosition="centertop",...
                Position=[10,235,220,320],...
                BackgroundColor=[0.90,0.90,0.90],Enable="off");

            obj.PolaritySwitch = uiswitch(obj.Container,...
                Items=["positive","negative"],...
                ItemsData=["positive","negative"],...
                Position=[85,270,45,20]);
            obj.PolaritySwitch.Value = CallingData.MSPolarity;
            obj.PolLabel  = uilabel(obj.Container,"Text","MS Polarity",...
                "Position",[75,250,70,20],...
                "HorizontalAlignment","center","VerticalAlignment","center");

            obj.BlankTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.BLKSubtraction,...
                "Position",[35,228,175,20],...
                "Text","Blank Subtraction");

            obj.PeakAlignmentTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.Peakalign,...
                "Position",[35,203,175,20],...
                "Text","Peak Alignment");
            

            obj.MSAlignmentTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.MSalign,...
                "Position",[35,178,175,20],...
                "Text","Spectral Alignment");
            

            obj.BaseCorrTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.BaseCorr,...
                "Position",[35,153,175,20],...
                "Text","Baseline Correction");
            

            obj.GolayTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.Smoothing,...
                "Position",[35,128,175,20],...
                "Text","Golay Smoothing");
            

            obj.ScalingTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.ScalingCorr,...
                "Position",[35,103,175,20],...
                "Text","Sample and Group Scaling");

            obj.ISTDTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.ISTDCorr,...
                "Position",[35,78,175,20],...
                "Text","ISTD Normalization");
            

            obj.IsotopeFilterTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.IsotopeFilter,...
                "Position",[35,53,175,20],...
                "Text","Isotope Filter");
            

            obj.AddNLossTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.IsotopeFilter,...S
                "Position",[35,28,175,20],...
                "Text","Adduct & Neutral Loss Filter");
            

            obj.ContaminantTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.ContaminantFilter,...
                "Position",[35,3,175,20],...
                "Text","Common Contaminant Filter");
            
           
            
        end
    end
        
end