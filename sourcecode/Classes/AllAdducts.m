classdef AllAdducts < AddNLoss
    %% Object for storing Negative Adduct Filter UI elements
    properties
        SingleArray matlab.ui.control.CheckBox
        AllSingles matlab.ui.control.CheckBox
        DoubleArray matlab.ui.control.CheckBox
        AllDoubles matlab.ui.control.CheckBox
        DimerArray matlab.ui.control.CheckBox
        AllDimers matlab.ui.control.CheckBox
        TripleArray matlab.ui.control.CheckBox
        AllTriples matlab.ui.control.CheckBox
    end

    methods
        function obj = AllAdducts(CallingData)

            obj@AddNLoss(CallingData)
            % Change Size of GUI to accomodate dual polarity adducts
            obj.UIFigure.Position = [100 100 815 600];
            obj.AdductsMainPanel.Position = [15 10 430 520];
            obj.SingleSubPanel.Position = [10 185 270 300];
            obj.DualSubPanel.Position = [285 185 135 300];
            obj.DimerSubPanel.Position = [10 5 270 170];
            obj.TripleSubPanel.Position = [285 5 135 170];
            obj.NLossMainPanel.Position = [460,70,350,460];
            obj.AcceptButton.Position =  [625 15 180 40];
            obj.CosSimField.Position = [550 30 50 20];
            obj.CosSimLabel.Position = [460 30 85 30];
            %Construct All Tickboxes
            %% Single Charged Adducts
            % Create S1
            obj.SingleArray(1) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(1).Text = 'M+NH4';
            obj.SingleArray(1).Position = [15 235 90 22];

            % Create S2
            obj.SingleArray(2) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(2).Text = 'M+Na';
            obj.SingleArray(2).Position = [15 215 90 22];

            % Create S3
            obj.SingleArray(3) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(3).Text = 'M+K';
            obj.SingleArray(3).Position = [15 195 90 22];

            % Create S4
            obj.SingleArray(4) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(4).Text = 'M+ACN+H';
            obj.SingleArray(4).Position = [15 175 90 22];

            % Create S5
            obj.SingleArray(5) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(5).Text = 'M+MeOH+H';
            obj.SingleArray(5).Position = [15 155 90 22];

            % Create S6
            obj.SingleArray(6) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(6).Text = 'M+IPA+H';
            obj.SingleArray(6).Position = [15 135 90 22];

            % Create S7
            obj.SingleArray(7) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(7).Text = 'M+DMSO+H';
            obj.SingleArray(7).Position = [15 115 90 22];

            % Create S8
            obj.SingleArray(8) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(8).Text = 'M+ACN+Na';
            obj.SingleArray(8).Position = [15 95 90 22];

            % Create S9
            obj.SingleArray(9) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(9).Text = 'M+2ACN+H';
            obj.SingleArray(9).Position = [15 75 90 22];

            % Create S10
            obj.SingleArray(10) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(10).Text = 'M+IPA+Na';
            obj.SingleArray(10).Position = [15 55 90 22];

            % Create S11
            obj.SingleArray(11) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(11).Text = 'M+2Na-H';
            obj.SingleArray(11).Position = [15 35 90 22];

            % Create S12
            obj.SingleArray(12) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(12).Text = 'M+2K-H';
            obj.SingleArray(12).Position = [15 15 90 22];

            % Create S13
            obj.SingleArray(13) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(13).Text = 'M+Cl';
            obj.SingleArray(13).Position = [150 235 90 22];

            % Create S14
            obj.SingleArray(14) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(14).Text = 'M+Br';
            obj.SingleArray(14).Position = [150 215 90 22];

            % Create S15
            obj.SingleArray(15) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(15).Text = 'M+Fac-H';
            obj.SingleArray(15).Position = [150 195 90 22];

            % Create S16
            obj.SingleArray(16) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(16).Text = 'M+HAc-H';
            obj.SingleArray(16).Position = [150 175 90 22];

            % Create S17
            obj.SingleArray(17) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(17).Text = 'M+NaAc-H';
            obj.SingleArray(17).Position = [150 155 90 22];

            % Create S18
            obj.SingleArray(18) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(18).Text = 'M+NH4Ac-H';
            obj.SingleArray(18).Position = [150 135 90 22];

            % Create S19
            obj.SingleArray(19) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(19).Text = 'M+TFA-H';
            obj.SingleArray(19).Position = [150 115 90 22];

            % Create S20
            obj.SingleArray(20) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(20).Text = 'M+Na-2H';
            obj.SingleArray(20).Position = [150 95 90 22];

            % Create S21
            obj.SingleArray(21) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(21).Text = 'M+Na-2H';
            obj.SingleArray(21).Position = [150 75 90 22];

            % Create S22
            obj.SingleArray(22) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(22).Text = 'M+K-2H';
            obj.SingleArray(22).Position = [150 55 90 22];

            % Create AllSingle
            obj.AllSingles = uicheckbox(obj.SingleSubPanel);
            obj.AllSingles.ValueChangedFcn = @(src,event) {AllSingleChanged(obj,src,event)};
            obj.AllSingles.Text = 'All';
            obj.AllSingles.Position = [15 255 54 22];

            %% Dimers
            % Create D1
            obj.DimerArray(1) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(1).Text = '2M+H';
            obj.DimerArray(1).Position = [15 105 90 22];

            % Create D2
            obj.DimerArray(2) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(2).Text = '2M+NH4';
            obj.DimerArray(2).Position = [15 85 90 22];

            % Create D3
            obj.DimerArray(3) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(3).Text = '2M+Na';
            obj.DimerArray(3).Position = [15 65 90 22];

            % Create D4
            obj.DimerArray(4) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(4).Text = '2M+K';
            obj.DimerArray(4).Position = [15 45 90 22];

            % Create D5
            obj.DimerArray(5) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(5).Text = '2M+ACN+H';
            obj.DimerArray(5).Position = [15 25 90 22];

            % Create D6
            obj.DimerArray(6) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(6).Text = '2M+ACN+Na';
            obj.DimerArray(6).Position = [15 5 90 22];

            % Create D7
            obj.DimerArray(7) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(7).Text = '2M-H';
            obj.DimerArray(7).Position = [150 105 90 22];

            % Create D8
            obj.DimerArray(8) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(8).Text = '2M+Fac-H';
            obj.DimerArray(8).Position = [150 85 90 22];

            % Create D9
            obj.DimerArray(9) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(9).Text = '2M+HAc-H';
            obj.DimerArray(9).Position = [150 65 90 22];

            % Create D10
            obj.DimerArray(10) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(10).Text = '3M-H';
            obj.DimerArray(10).Position = [150 45 90 22];

            % Create AllD
            obj.AllDimers = uicheckbox(obj.DimerSubPanel);
            obj.AllDimers.ValueChangedFcn = @(src,event) {AllDimerChanged(obj,src,event)};
            obj.AllDimers.Text = 'All';
            obj.AllDimers.Position = [15 125 90 22];

            %% Triple Charged
            % Create T1
            obj.TripleArray(1) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(1).Text = 'M+3H';
            obj.TripleArray(1).Position = [15 105 90 22];

            % Create T2
            obj.TripleArray(2) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(2).Text = 'M+Na+2H';
            obj.TripleArray(2).Position = [15 85 90 22];

            % Create T3
            obj.TripleArray(3) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(3).Text = 'M+2Na+H';
            obj.TripleArray(3).Position = [15 65 90 22];

            % Create T4
            obj.TripleArray(4) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(4).Text = 'M+3Na';
            obj.TripleArray(4).Position = [15 45 90 22];

            % Create T5
            obj.TripleArray(5) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(5).Text = 'M-3H';
            obj.TripleArray(5).Position = [15 25 90 22];

            % Create All Triples
            obj.AllTriples = uicheckbox(obj.TripleSubPanel);
            obj.AllTriples.ValueChangedFcn = @(src,event) {AllTriplesChanged(obj,src,event)};
            obj.AllTriples.Text = 'All';
            obj.AllTriples.Position = [15 125 35 22];

            %% Dual Charged Adducts
            % Create Dual1
            obj.DoubleArray(1) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(1).Text = 'M+2H';
            obj.DoubleArray(1).Position = [15 235 54 22];

            % Create Dual2
            obj.DoubleArray(2) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(2).Text = 'M+NH4+H';
            obj.DoubleArray(2).Position = [15 215 79 22];

            % Create Dual3
            obj.DoubleArray(3) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(3).Text = 'M+Na+H';
            obj.DoubleArray(3).Position = [15 195 70 22];

            % Create Dual4
            obj.DoubleArray(4) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(4).Text = 'M+K+H';
            obj.DoubleArray(4).Position = [15 175 63 22];

            % Create Dual5
            obj.DoubleArray(5) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(5).Text = 'M+2Na';
            obj.DoubleArray(5).Position = [15 155 61 22];

            % Create Dual6
            obj.DoubleArray(6) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(6).Text = 'M+ACN+2H';
            obj.DoubleArray(6).Position = [15 135 87 22];

            % Create Dual7
            obj.DoubleArray(7) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(7).Text = 'M+2ACN+2H';
            obj.DoubleArray(7).Position = [15 115 93 22];

            % Create Dual8
            obj.DoubleArray(8) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(8).Text = 'M+3ACN+2H';
            obj.DoubleArray(8).Position = [15 95 93 22];

            % Create Dual9
            obj.DoubleArray(9) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(9).Text = 'M-2H';
            obj.DoubleArray(9).Position = [15 75 93 22];

            % Create All Dual
            obj.AllDoubles = uicheckbox(obj.DualSubPanel);
            obj.AllDoubles.ValueChangedFcn = @(src,event) {AllDualsChanged(obj,src,event)};
            obj.AllDoubles.Text = 'All';
            obj.AllDoubles.Position = [15 255 35 22];

            %load previous pos data
            val = num2cell(CallingData.AddSelectedPos(27:30));
            [obj.TripleArray(1:4).Value] = deal(val{:});
            val = num2cell(CallingData.AddSelectedPos(19:26));
            [obj.DoubleArray(1:8).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedPos(13:18));
            [obj.DimerArray(1:6).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedPos(1:12));
            [obj.SingleArray(1:12).Value] = deal(val{:});

            %load previous neg data
            obj.TripleArray(5).Value = CallingData.AddSelectedNeg(16);
            obj.DoubleArray(9).Value = CallingData.AddSelectedNeg(15);
            val=num2cell(CallingData.AddSelectedNeg(11:14));
            [obj.DimerArray(7:10).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedNeg(1:10));
            [obj.SingleArray(13:22).Value] = deal(val{:});
        end

        function obj = AllSingleChanged(obj,~,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.SingleArray,2));
            [obj.SingleArray(:).Value] = val{:};
        end
        function obj = AllDimerChanged(obj,~,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.DimerArray,2));
            [obj.DimerArray(:).Value] = val{:};
        end
        function obj = AllTriplesChanged(obj,~,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.TripleArray,2));
            [obj.TripleArray(:).Value] = val{:};
        end
        function obj = AllDualsChanged(obj,~,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.DoubleArray,2));
            [obj.DoubleArray(:).Value] = val{:};
        end
    end
end