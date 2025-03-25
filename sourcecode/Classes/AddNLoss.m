classdef AddNLoss
    %% Superclass for Adduct filter objects
    properties
        uiFigure matlab.ui.Figure
        label
        cosSimField
        cosSimLabel
        adductsMainPanel
        singleSubPanel
        dualSubPanel
        dimerSubPanel
        tripleSubPanel
        nLossMainPanel
        smolSubPanel
        conjugateSubPanel
        acceptButton
        smolArray matlab.ui.control.CheckBox
        conjugatesArray matlab.ui.control.CheckBox
        allSmol matlab.ui.control.CheckBox
        allConjungates matlab.ui.control.CheckBox
    end

    methods
        function obj = AddNLoss(CallingData)

            if nargin > 0
                obj.uiFigure = uifigure(Position=[100,100,680,600],Name="Adduct and Neutral Loss Selection",Visible="off");
                obj.label = uilabel(obj.uiFigure,"Text","Define Adducts and Neutral losses",...
                    "Position",[15,545,540,45],...
                    FontSize=32,...
                    FontWeight="bold");
                % Construct Adduct Main and Sub Panels
                obj.adductsMainPanel = uipanel(obj.uiFigure,"Title","Adducts",...
                    "BackgroundColor",[0.9,0.9,0.9],...
                    "Position",[15,10,300,520],...
                    FontWeight="bold",...
                    FontSize=18);
                obj.singleSubPanel = uipanel(obj.adductsMainPanel,"Title","Single Charged", ...
                    "Position",[10,185,135,300]);
                obj.dualSubPanel = uipanel(obj.adductsMainPanel,"Title","Dual Charged", ...
                    "Position",[150,185,135,300]);
                obj.dimerSubPanel = uipanel(obj.adductsMainPanel,"Title","Analyte Dimer", ...
                    "Position",[10,5,135,170]);
                obj.tripleSubPanel = uipanel(obj.adductsMainPanel,"Title","Triple Charged", ...
                    "Position",[150,5,135,170]);

                % Construct Neutral Loss Main and Sub Panels
                obj.nLossMainPanel = uipanel(obj.uiFigure,"Title","Neutral Loss",...
                    "BackgroundColor",[0.9,0.9,0.9],...
                    "Position",[325,70,350,460],...
                    FontWeight="bold",...
                    FontSize=18);
                obj.smolSubPanel = uipanel(obj.nLossMainPanel,"Title","Small Molecule",...
                    "Position",[10,10,110,410]);
                obj.conjugateSubPanel = uipanel(obj.nLossMainPanel,"Title","Conjungates",...
                    "Position",[125,10,210,410]);

                % Construct Cosine Similarity field and Button
                obj.acceptButton = uibutton(obj.uiFigure,...
                    "Text","Accept", ...
                    "Position",[490,15,180,40],...
                    "FontSize",22,...
                    FontWeight="bold");

                obj.cosSimField = uieditfield(obj.uiFigure,"numeric",...
                    "ValueDisplayFormat","%.1f %%",...
                    Position=[415,30,50,20],...
                    Limits=[0,100],...
                    Value=CallingData.CosSim*100,...
                    Tooltip="Minimum peak shape similarity between original peak and possible adduct peak. Calculated as cosine similarity");
                obj.cosSimLabel = uilabel(obj.uiFigure,"Text",["minimum","Peak Similarity"],...
                    "Position",[325,30,85,30],...
                    "HorizontalAlignment","left");

                % Construct Common TickBoxes
                % Small Molecule
                % Create Smol1
                obj.smolArray(1) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(1).Text = '-NH3';
                obj.smolArray(1).Position = [15 345 50 22];

                % Create Smol2
                obj.smolArray(2) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(2).Text = '-H2O';
                obj.smolArray(2).Position = [15 325 51 22];

                % Create Smol3
                obj.smolArray(3) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(3).Text = '-HCN';
                obj.smolArray(3).Position = [15 305 52 22];

                % Create Smol4
                obj.smolArray(4) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(4).Text = '-CO';
                obj.smolArray(4).Position = [15 285 44 22];

                % Create Smol5
                obj.smolArray(5) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(5).Text = '-H2CN';
                obj.smolArray(5).Position = [15 265 59 22];

                % Create Smol6
                obj.smolArray(6) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(6).Text = '-NO';
                obj.smolArray(6).Position = [15 245 44 22];

                % Create Smol7
                obj.smolArray(7) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(7).Text = '-CH2O';
                obj.smolArray(7).Position = [15 225 59 22];

                % Create Smol8
                obj.smolArray(8) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(8).Text = '-CH3OH';
                obj.smolArray(8).Position = [15 205 68 22];

                % Create Smol9
                obj.smolArray(9) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(9).Text = '-H2S';
                obj.smolArray(9).Position = [15 185 49 22];

                % Create Smol10
                obj.smolArray(10) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(10).Text = '-HCl';
                obj.smolArray(10).Position = [15 165 46 22];

                % Create Smol11
                obj.smolArray(11) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(11).Text = '-CO2';
                obj.smolArray(11).Position = [15 145 51 22];

                % Create Smol12
                obj.smolArray(12) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(12).Text = '-NO2';
                obj.smolArray(12).Position = [15 125 51 22];

                % Create Smol13
                obj.smolArray(13) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(13).Text = '-HCOOH';
                obj.smolArray(13).Position = [15 105 71 22];

                % Create Smol14
                obj.smolArray(14) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(14).Text = '-SO2';
                obj.smolArray(14).Position = [15 85 50 22];

                % Create Smol15
                obj.smolArray(15) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(15).Text = '-HPO3';
                obj.smolArray(15).Position = [15 65 59 22];

                % Create Smol16
                obj.smolArray(16) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(16).Text = '-H2SO3';
                obj.smolArray(16).Position = [15 45 65 22];

                % Create Smol17
                obj.smolArray(17) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(17).Text = '-H2PO4';
                obj.smolArray(17).Position = [15 25 65 22];

                % Create Smol18
                obj.smolArray(18) = uicheckbox(obj.smolSubPanel);
                obj.smolArray(18).Text = '-HI';
                obj.smolArray(18).Position = [15 5 38 22];

                % Create AllSmol
                obj.allSmol = uicheckbox(obj.smolSubPanel);
                obj.allSmol.ValueChangedFcn = @(src,event) {allSmolChanged(obj,src,event)};
                obj.allSmol.Text = 'All';
                obj.allSmol.Position = [15 365 35 22];

                % Conjungates
                % Create C1
                obj.conjugatesArray(1) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(1).Text = '-Histidine';
                obj.conjugatesArray(1).Position = [15 345 72 22];

                % Create C2
                obj.conjugatesArray(2) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(2).Text = '-Cysteine';
                obj.conjugatesArray(2).Position = [15 325 73 22];

                % Create C3
                obj.conjugatesArray(3) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(3).Text = '-Didesoxyhexoside';
                obj.conjugatesArray(3).Position = [15 305 123 22];

                % Create C4
                obj.conjugatesArray(4) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(4).Text = '-Pentoside';
                obj.conjugatesArray(4).Position = [15 285 79 22];

                % Create C5
                obj.conjugatesArray(5) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(5).Text = '-Deoxyhexoside';
                obj.conjugatesArray(5).Position = [15 265 108 22];

                % Create C6
                obj.conjugatesArray(6) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(6).Text = '-Glutamic Acid';
                obj.conjugatesArray(6).Position = [15 245 99 22];

                % Create C7
                obj.conjugatesArray(7) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(7).Text = '-Hexoside';
                obj.conjugatesArray(7).Position = [15 225 76 22];

                % Create C8
                obj.conjugatesArray(8) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(8).Text = '-N-acetylcysteine';
                obj.conjugatesArray(8).Position = [15 205 114 22];

                % Create C9
                obj.conjugatesArray(9) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(9).Text = '-Rhamnoside';
                obj.conjugatesArray(9).Position = [15 185 93 22];

                % Create C10
                obj.conjugatesArray(10) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(10).Text = '-Glucuronide';
                obj.conjugatesArray(10).Position = [15 165 91 22];

                % Create C11
                obj.conjugatesArray(11) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(11).Text = '-Glucuronide (benzylic)';
                obj.conjugatesArray(11).Position = [15 145 145 22];

                % Create C12
                obj.conjugatesArray(12) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(12).Text = '-N-acetylglucosamine (benzylic)';
                obj.conjugatesArray(12).Position = [15 125 193 22];

                % Create C13
                obj.conjugatesArray(13) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(13).Text = '-N-acetylglucosamine';
                obj.conjugatesArray(13).Position = [15 105 137 22];

                % Create C14
                obj.conjugatesArray(14) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(14).Text = '-Malonylglucuronide';
                obj.conjugatesArray(14).Position = [15 85 129 22];

                % Create C15
                obj.conjugatesArray(15) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(15).Text = '-gamma GluCys';
                obj.conjugatesArray(15).Position = [15 65 109 22];

                % Create C16
                obj.conjugatesArray(16) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(16).Text = '-Malonylglucuronide (benzylic)';
                obj.conjugatesArray(16).Position = [15 45 184 22];

                % Create C17
                obj.conjugatesArray(17) = uicheckbox(obj.conjugateSubPanel);
                obj.conjugatesArray(17).Text = '-Gluthatione';
                obj.conjugatesArray(17).Position = [15 25 87 22];

                % Create AllC
                obj.allConjungates = uicheckbox(obj.conjugateSubPanel);
                obj.allConjungates.ValueChangedFcn = @(src,event) {allConjungateChanged(obj,src,event)};
                obj.allConjungates.Text = 'All';
                obj.allConjungates.Position = [15 365 35 22];

                %load previous values
                val=num2cell(CallingData.NeutralSelectedCon);
                [obj.conjugatesArray(:).Value] = deal(val{:});
                val = num2cell(CallingData.NeutralSelectedSmol);
                [obj.smolArray(:).Value] = deal(val{:});
            end
        end

        function obj = allSmolChanged(obj,src,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.smolArray,2));
            [obj.smolArray(:).Value] = val{:};
        end
        function obj = allConjungateChanged(obj,src,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.conjugatesArray,2));
            [obj.conjugatesArray(:).Value] = val{:};
        end
    end
end