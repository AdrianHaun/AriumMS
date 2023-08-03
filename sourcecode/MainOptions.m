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
                BackgroundColor=[0.90,0.90,0.90],Enable="off",...
                Tooltip="Main settings. Activating an option activates the corresponding panel.");

            obj.PolaritySwitch = uiswitch(obj.Container,...
                Items=["positive","negative"],...
                ItemsData=["positive","negative"],...
                Position=[85,270,45,20]);
            obj.PolaritySwitch.Value = CallingData.MSPolarity;
            obj.PolLabel  = uilabel(obj.Container,"Text","MS Polarity",...
                "Position",[75,250,70,20],...
                "HorizontalAlignment","center","VerticalAlignment","center",...
                Tooltip="Set the polarity of the measurement. This setting affects the mass list used for the adduct and contaminant filter");

            obj.BlankTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.BLKSubtraction,...
                "Position",[35,228,175,20],...
                "Text","Blank Subtraction",...
                Tooltip="Activates the subtraction of blank files. If several blanks are selected, an average is generated first.");

            obj.PeakAlignmentTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.Peakalign,...
                "Position",[35,203,175,20],...
                "Text","Peak Alignment",...
                Tooltip="Aligns peaks between files to reference peaks. For each ROI the highest signal is used as reference peak. First, it creates a synthetic signal from the reference peaks using Gaussian pulses centered at the separation-unit values specified by the reference peak. Then, it shifts and scales the separation-unit scale to find the maximum alignment between the input signals and the synthetic signal. (It uses an iterative multiresolution grid search until it finds the best scale and shift factors for each signal.) Once the new separation-unit scale is determined, the corrected signals are created by resampling their intensities at the original separation-unit values. The resampling method preserves the shape of the peaks.");
            

            obj.MSAlignmentTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.MSalign,...
                "Position",[35,178,175,20],...
                "Text","Spectral Alignment",...
                Tooltip="Aligns mass spectra from multiple peak lists (centroided data) by first estimating a list of common m/z (CMZ vector) values estimated by considering the peaks in all spectra in a file. It then aligns the peaks in each spectrum to the CMZ vector. Only recommended for low resolution measurements");
            

            obj.BaseCorrTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.BaseCorr,...
                "Position",[35,153,175,20],...
                "Text","Baseline Correction",...
                Tooltip="Corrects variable baseline by estimating the baseline within multiple shifted windows, regressing the varying baseline to the window points using a spline approximation, then adjusting the baseline of the peak signals.");
            

            obj.GolayTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.Smoothing,...
                "Position",[35,128,175,20],...
                "Text","Golay Smoothing",...
                Tooltip="Smoothing of the signals using a Savitzky-Golay filter.");
            

            obj.ScalingTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.ScalingCorr,...
                "Position",[35,103,175,20],...
                "Text","Sample and Group Scaling",...
                Tooltip="Scale the feature intensities of this group. Specify a scaling factor for the whole group and/or scaling factors for each individual measurement.");

            obj.ISTDTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.ISTDCorr,...
                "Position",[35,78,175,20],...
                "Text","ISTD Normalization",...
                Tooltip="Normalize feature intensities and correct mass axis with internal standard(s)");
            

            obj.IsotopeFilterTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.IsotopeFilter,...
                "Position",[35,53,175,20],...
                "Text","Isotope Filter",...
                Tooltip="Removes isotopic peaks. Verified are peak shape, intensity related to isotope abundance, retention time and mass tolerance.");
            

            obj.AddNLossTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.IsotopeFilter,...S
                "Position",[35,28,175,20],...
                "Text","Adduct & Neutral Loss Filter",...
                Tooltip="Removes adduct and neutral loss peaks. Verified are peak shape, retention time and mass tolerance.");
            

            obj.ContaminantTickBox = uicheckbox(obj.Container,...
                "Value",CallingData.ContaminantFilter,...
                "Position",[35,3,175,20],...
                "Text","Common Contaminant Filter",...
                Tooltip="Removes all features with masses listed in the UWPR Common Mass Spec Contaminants list.");
            
           
            
        end
    end
        
end