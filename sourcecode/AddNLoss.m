classdef AddNLoss
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        UIFigure matlab.ui.Figure
        Label
        CosSimField
        CosSimLabel
        AdductsMainPanel
        SingleSubPanel
        DualSubPanel
        DimerSubPanel
        TripleSubPanel
        NLossMainPanel
        SmolSubPanel
        ConjugateSubPanel
        AcceptButton
        SmolArray matlab.ui.control.CheckBox
        ConjugatesArray matlab.ui.control.CheckBox
        AllSmol matlab.ui.control.CheckBox
        AllConjungates matlab.ui.control.CheckBox
    end

    methods
        function obj = AddNLoss(CallingData)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here

            obj.UIFigure = uifigure(Position=[100,100,680,600],Name="Adduct and Neutral Loss Selection",Visible="off");
            obj.Label = uilabel(obj.UIFigure,"Text","Define Adducts and Neutral losses",...
                "Position",[15,545,540,45],...
                FontSize=32,...
                FontWeight="bold");
            % Construct Adduct Main and Sub Panels
            obj.AdductsMainPanel = uipanel(obj.UIFigure,"Title","Adducts",...
                "BackgroundColor",[0.9,0.9,0.9],...
                "Position",[15,10,300,520],...
                FontWeight="bold",...
                FontSize=18);
            obj.SingleSubPanel = uipanel(obj.AdductsMainPanel,"Title","Single Charged", ...
                "Position",[10,185,135,300]);
            obj.DualSubPanel = uipanel(obj.AdductsMainPanel,"Title","Dual Charged", ...
                "Position",[150,185,135,300]);
            obj.DimerSubPanel = uipanel(obj.AdductsMainPanel,"Title","Analyte Dimer", ...
                "Position",[10,5,135,170]);
            obj.TripleSubPanel = uipanel(obj.AdductsMainPanel,"Title","Triple Charged", ...
                "Position",[150,5,135,170]);

            % Construct Neutral Loss Main and Sub Panels
            obj.NLossMainPanel = uipanel(obj.UIFigure,"Title","Neutral Loss",...
                "BackgroundColor",[0.9,0.9,0.9],...
                "Position",[325,70,350,460],...
                FontWeight="bold",...
                FontSize=18);
            obj.SmolSubPanel = uipanel(obj.NLossMainPanel,"Title","Small Molecule",...
                "Position",[10,10,110,410]);
            obj.ConjugateSubPanel = uipanel(obj.NLossMainPanel,"Title","Conjungates",...
                "Position",[125,10,210,410]);

            % Construct Cosine Similarity field and Button
            obj.AcceptButton = uibutton(obj.UIFigure,...
                "Text","Accept", ...
                "Position",[490,15,180,40],...
                "FontSize",22,...
                FontWeight="bold");

            obj.CosSimField = uieditfield(obj.UIFigure,"numeric",...
                "ValueDisplayFormat","%.1f %%",...
                Position=[415,30,50,20],...
                Limits=[0,100],...
                Value=CallingData.CosSim*100,...
                Tooltip="Minimum peak shape similarity between original peak and possible adduct peak. Calculated as cosine similarity");
            obj.CosSimLabel = uilabel(obj.UIFigure,"Text",["minimum","Peak Similarity"],...
                "Position",[325,30,85,30],...
                "HorizontalAlignment","left");

            % Construct Common TickBoxes
            % Small Molecule
            % Create Smol1
            obj.SmolArray(1) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(1).Text = '-NH3';
            obj.SmolArray(1).Position = [15 345 50 22];

            % Create Smol2
            obj.SmolArray(2) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(2).Text = '-H2O';
            obj.SmolArray(2).Position = [15 325 51 22];

            % Create Smol3
            obj.SmolArray(3) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(3).Text = '-HCN';
            obj.SmolArray(3).Position = [15 305 52 22];

            % Create Smol4
            obj.SmolArray(4) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(4).Text = '-CO';
            obj.SmolArray(4).Position = [15 285 44 22];

            % Create Smol5
            obj.SmolArray(5) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(5).Text = '-H2CN';
            obj.SmolArray(5).Position = [15 265 59 22];

            % Create Smol6
            obj.SmolArray(6) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(6).Text = '-NO';
            obj.SmolArray(6).Position = [15 245 44 22];

            % Create Smol7
            obj.SmolArray(7) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(7).Text = '-CH2O';
            obj.SmolArray(7).Position = [15 225 59 22];

            % Create Smol8
            obj.SmolArray(8) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(8).Text = '-CH3OH';
            obj.SmolArray(8).Position = [15 205 68 22];

            % Create Smol9
            obj.SmolArray(9) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(9).Text = '-H2S';
            obj.SmolArray(9).Position = [15 185 49 22];

            % Create Smol10
            obj.SmolArray(10) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(10).Text = '-HCl';
            obj.SmolArray(10).Position = [15 165 46 22];

            % Create Smol11
            obj.SmolArray(11) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(11).Text = '-CO2';
            obj.SmolArray(11).Position = [15 145 51 22];

            % Create Smol12
            obj.SmolArray(12) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(12).Text = '-NO2';
            obj.SmolArray(12).Position = [15 125 51 22];

            % Create Smol13
            obj.SmolArray(13) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(13).Text = '-HCOOH';
            obj.SmolArray(13).Position = [15 105 71 22];

            % Create Smol14
            obj.SmolArray(14) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(14).Text = '-SO2';
            obj.SmolArray(14).Position = [15 85 50 22];

            % Create Smol15
            obj.SmolArray(15) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(15).Text = '-HPO3';
            obj.SmolArray(15).Position = [15 65 59 22];

            % Create Smol16
            obj.SmolArray(16) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(16).Text = '-H2SO3';
            obj.SmolArray(16).Position = [15 45 65 22];

            % Create Smol17
            obj.SmolArray(17) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(17).Text = '-H2PO4';
            obj.SmolArray(17).Position = [15 25 65 22];

            % Create Smol18
            obj.SmolArray(18) = uicheckbox(obj.SmolSubPanel);
            obj.SmolArray(18).Text = '-HI';
            obj.SmolArray(18).Position = [15 5 38 22];

            % Create AllSmol
            obj.AllSmol = uicheckbox(obj.SmolSubPanel);
            obj.AllSmol.ValueChangedFcn = @(src,event) {AllSmolChanged(obj,src,event)};
            obj.AllSmol.Text = 'All';
            obj.AllSmol.Position = [15 365 35 22];

            % Conjungates
            % Create C1
            obj.ConjugatesArray(1) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(1).Text = '-Histidine';
            obj.ConjugatesArray(1).Position = [15 345 72 22];

            % Create C2
            obj.ConjugatesArray(2) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(2).Text = '-Cysteine';
            obj.ConjugatesArray(2).Position = [15 325 73 22];

            % Create C3
            obj.ConjugatesArray(3) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(3).Text = '-Didesoxyhexoside';
            obj.ConjugatesArray(3).Position = [15 305 123 22];

            % Create C4
            obj.ConjugatesArray(4) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(4).Text = '-Pentoside';
            obj.ConjugatesArray(4).Position = [15 285 79 22];

            % Create C5
            obj.ConjugatesArray(5) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(5).Text = '-Deoxyhexoside';
            obj.ConjugatesArray(5).Position = [15 265 108 22];

            % Create C6
            obj.ConjugatesArray(6) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(6).Text = '-Glutamic Acid';
            obj.ConjugatesArray(6).Position = [15 245 99 22];

            % Create C7
            obj.ConjugatesArray(7) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(7).Text = '-Hexoside';
            obj.ConjugatesArray(7).Position = [15 225 76 22];

            % Create C8
            obj.ConjugatesArray(8) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(8).Text = '-N-acetylcysteine';
            obj.ConjugatesArray(8).Position = [15 205 114 22];

            % Create C9
            obj.ConjugatesArray(9) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(9).Text = '-Rhamnoside';
            obj.ConjugatesArray(9).Position = [15 185 93 22];

            % Create C10
            obj.ConjugatesArray(10) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(10).Text = '-Glucuronide';
            obj.ConjugatesArray(10).Position = [15 165 91 22];

            % Create C11
            obj.ConjugatesArray(11) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(11).Text = '-Glucuronide (benzylic)';
            obj.ConjugatesArray(11).Position = [15 145 145 22];

            % Create C12
            obj.ConjugatesArray(12) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(12).Text = '-N-acetylglucosamine (benzylic)';
            obj.ConjugatesArray(12).Position = [15 125 193 22];

            % Create C13
            obj.ConjugatesArray(13) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(13).Text = '-N-acetylglucosamine';
            obj.ConjugatesArray(13).Position = [15 105 137 22];

            % Create C14
            obj.ConjugatesArray(14) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(14).Text = '-Malonylglucuronide';
            obj.ConjugatesArray(14).Position = [15 85 129 22];

            % Create C15
            obj.ConjugatesArray(15) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(15).Text = '-gamma GluCys';
            obj.ConjugatesArray(15).Position = [15 65 109 22];

            % Create C16
            obj.ConjugatesArray(16) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(16).Text = '-Malonylglucuronide (benzylic)';
            obj.ConjugatesArray(16).Position = [15 45 184 22];

            % Create C17
            obj.ConjugatesArray(17) = uicheckbox(obj.ConjugateSubPanel);
            obj.ConjugatesArray(17).Text = '-Gluthatione';
            obj.ConjugatesArray(17).Position = [15 25 87 22];

            % Create AllC
            obj.AllConjungates = uicheckbox(obj.ConjugateSubPanel);
            obj.AllConjungates.ValueChangedFcn = @(src,event) {AllCChanged(obj,src,event)};
            obj.AllConjungates.Text = 'All';
            obj.AllConjungates.Position = [15 365 35 22];

            %load previous values
            val=num2cell(CallingData.NeutralSelectedCon);
            [obj.ConjugatesArray(:).Value] = deal(val{:});
            val = num2cell(CallingData.NeutralSelectedSmol);
            [obj.SmolArray(:).Value] = deal(val{:});
        end

        function obj = AllSmolChanged(obj,src,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.SmolArray,2));
            [obj.SmolArray(:).Value] = val{:};
        end
        function obj = AllCChanged(obj,src,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.ConjugatesArray,2));
            [obj.ConjugatesArray(:).Value] = val{:};
        end
    end
end