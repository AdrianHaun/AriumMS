classdef AllAdducts < AddNLoss
    %% Object for storing Negative Adduct Filter UI elements
    properties
        singleArray matlab.ui.control.CheckBox
        allSingles matlab.ui.control.CheckBox
        doubleArray matlab.ui.control.CheckBox
        allDoubles matlab.ui.control.CheckBox
        dimerArray matlab.ui.control.CheckBox
        allDimers matlab.ui.control.CheckBox
        tripleArray matlab.ui.control.CheckBox
        allTriples matlab.ui.control.CheckBox
    end

    methods
        function obj = AllAdducts(CallingData)

            obj@AddNLoss(CallingData)
            % Change Size of GUI to accommodate dual polarity adducts
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
            %Construct All Tick boxes
            %% Single Charged Adducts
            % Create S1
            obj.singleArray(1) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(1).Text = 'M+NH4';
            obj.singleArray(1).Position = [15 235 90 22];

            % Create S2
            obj.singleArray(2) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(2).Text = 'M+Na';
            obj.singleArray(2).Position = [15 215 90 22];

            % Create S3
            obj.singleArray(3) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(3).Text = 'M+K';
            obj.singleArray(3).Position = [15 195 90 22];

            % Create S4
            obj.singleArray(4) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(4).Text = 'M+ACN+H';
            obj.singleArray(4).Position = [15 175 90 22];

            % Create S5
            obj.singleArray(5) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(5).Text = 'M+MeOH+H';
            obj.singleArray(5).Position = [15 155 90 22];

            % Create S6
            obj.singleArray(6) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(6).Text = 'M+IPA+H';
            obj.singleArray(6).Position = [15 135 90 22];

            % Create S7
            obj.singleArray(7) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(7).Text = 'M+DMSO+H';
            obj.singleArray(7).Position = [15 115 90 22];

            % Create S8
            obj.singleArray(8) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(8).Text = 'M+ACN+Na';
            obj.singleArray(8).Position = [15 95 90 22];

            % Create S9
            obj.singleArray(9) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(9).Text = 'M+2ACN+H';
            obj.singleArray(9).Position = [15 75 90 22];

            % Create S10
            obj.singleArray(10) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(10).Text = 'M+IPA+Na';
            obj.singleArray(10).Position = [15 55 90 22];

            % Create S11
            obj.singleArray(11) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(11).Text = 'M+2Na-H';
            obj.singleArray(11).Position = [15 35 90 22];

            % Create S12
            obj.singleArray(12) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(12).Text = 'M+2K-H';
            obj.singleArray(12).Position = [15 15 90 22];

            % Create S13
            obj.singleArray(13) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(13).Text = 'M+Cl';
            obj.singleArray(13).Position = [150 235 90 22];

            % Create S14
            obj.singleArray(14) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(14).Text = 'M+Br';
            obj.singleArray(14).Position = [150 215 90 22];

            % Create S15
            obj.singleArray(15) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(15).Text = 'M+Fac-H';
            obj.singleArray(15).Position = [150 195 90 22];

            % Create S16
            obj.singleArray(16) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(16).Text = 'M+HAc-H';
            obj.singleArray(16).Position = [150 175 90 22];

            % Create S17
            obj.singleArray(17) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(17).Text = 'M+NaAc-H';
            obj.singleArray(17).Position = [150 155 90 22];

            % Create S18
            obj.singleArray(18) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(18).Text = 'M+NH4Ac-H';
            obj.singleArray(18).Position = [150 135 90 22];

            % Create S19
            obj.singleArray(19) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(19).Text = 'M+TFA-H';
            obj.singleArray(19).Position = [150 115 90 22];

            % Create S20
            obj.singleArray(20) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(20).Text = 'M+Na-2H';
            obj.singleArray(20).Position = [150 95 90 22];

            % Create S21
            obj.singleArray(21) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(21).Text = 'M+Na-2H';
            obj.singleArray(21).Position = [150 75 90 22];

            % Create S22
            obj.singleArray(22) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(22).Text = 'M+K-2H';
            obj.singleArray(22).Position = [150 55 90 22];

            % Create AllSingle
            obj.allSingles = uicheckbox(obj.SingleSubPanel);
            obj.allSingles.ValueChangedFcn = @(src,event) {allSingleChanged(obj,src,event)};
            obj.allSingles.Text = 'All';
            obj.allSingles.Position = [15 255 54 22];

            %% Dimers
            % Create D1
            obj.dimerArray(1) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(1).Text = '2M+H';
            obj.dimerArray(1).Position = [15 105 90 22];

            % Create D2
            obj.dimerArray(2) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(2).Text = '2M+NH4';
            obj.dimerArray(2).Position = [15 85 90 22];

            % Create D3
            obj.dimerArray(3) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(3).Text = '2M+Na';
            obj.dimerArray(3).Position = [15 65 90 22];

            % Create D4
            obj.dimerArray(4) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(4).Text = '2M+K';
            obj.dimerArray(4).Position = [15 45 90 22];

            % Create D5
            obj.dimerArray(5) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(5).Text = '2M+ACN+H';
            obj.dimerArray(5).Position = [15 25 90 22];

            % Create D6
            obj.dimerArray(6) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(6).Text = '2M+ACN+Na';
            obj.dimerArray(6).Position = [15 5 90 22];

            % Create D7
            obj.dimerArray(7) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(7).Text = '2M-H';
            obj.dimerArray(7).Position = [150 105 90 22];

            % Create D8
            obj.dimerArray(8) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(8).Text = '2M+Fac-H';
            obj.dimerArray(8).Position = [150 85 90 22];

            % Create D9
            obj.dimerArray(9) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(9).Text = '2M+HAc-H';
            obj.dimerArray(9).Position = [150 65 90 22];

            % Create D10
            obj.dimerArray(10) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(10).Text = '3M-H';
            obj.dimerArray(10).Position = [150 45 90 22];

            % Create AllD
            obj.allDimers = uicheckbox(obj.DimerSubPanel);
            obj.allDimers.ValueChangedFcn = @(src,event) {allDimerChanged(obj,src,event)};
            obj.allDimers.Text = 'All';
            obj.allDimers.Position = [15 125 90 22];

            %% Triple Charged
            % Create T1
            obj.tripleArray(1) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(1).Text = 'M+3H';
            obj.tripleArray(1).Position = [15 105 90 22];

            % Create T2
            obj.tripleArray(2) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(2).Text = 'M+Na+2H';
            obj.tripleArray(2).Position = [15 85 90 22];

            % Create T3
            obj.tripleArray(3) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(3).Text = 'M+2Na+H';
            obj.tripleArray(3).Position = [15 65 90 22];

            % Create T4
            obj.tripleArray(4) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(4).Text = 'M+3Na';
            obj.tripleArray(4).Position = [15 45 90 22];

            % Create T5
            obj.tripleArray(5) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(5).Text = 'M-3H';
            obj.tripleArray(5).Position = [15 25 90 22];

            % Create All Triples
            obj.allTriples = uicheckbox(obj.TripleSubPanel);
            obj.allTriples.ValueChangedFcn = @(src,event) {allTriplesChanged(obj,src,event)};
            obj.allTriples.Text = 'All';
            obj.allTriples.Position = [15 125 35 22];

            %% Dual Charged Adducts
            % Create Dual1
            obj.doubleArray(1) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(1).Text = 'M+2H';
            obj.doubleArray(1).Position = [15 235 54 22];

            % Create Dual2
            obj.doubleArray(2) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(2).Text = 'M+NH4+H';
            obj.doubleArray(2).Position = [15 215 79 22];

            % Create Dual3
            obj.doubleArray(3) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(3).Text = 'M+Na+H';
            obj.doubleArray(3).Position = [15 195 70 22];

            % Create Dual4
            obj.doubleArray(4) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(4).Text = 'M+K+H';
            obj.doubleArray(4).Position = [15 175 63 22];

            % Create Dual5
            obj.doubleArray(5) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(5).Text = 'M+2Na';
            obj.doubleArray(5).Position = [15 155 61 22];

            % Create Dual6
            obj.doubleArray(6) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(6).Text = 'M+ACN+2H';
            obj.doubleArray(6).Position = [15 135 87 22];

            % Create Dual7
            obj.doubleArray(7) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(7).Text = 'M+2ACN+2H';
            obj.doubleArray(7).Position = [15 115 93 22];

            % Create Dual8
            obj.doubleArray(8) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(8).Text = 'M+3ACN+2H';
            obj.doubleArray(8).Position = [15 95 93 22];

            % Create Dual9
            obj.doubleArray(9) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(9).Text = 'M-2H';
            obj.doubleArray(9).Position = [15 75 93 22];

            % Create All Dual
            obj.allDoubles = uicheckbox(obj.DualSubPanel);
            obj.allDoubles.ValueChangedFcn = @(src,event) {allDualsChanged(obj,src,event)};
            obj.allDoubles.Text = 'All';
            obj.allDoubles.Position = [15 255 35 22];

            %load previous positive data
            val = num2cell(CallingData.AddSelectedPos(27:30));
            [obj.tripleArray(1:4).Value] = deal(val{:});
            val = num2cell(CallingData.AddSelectedPos(19:26));
            [obj.doubleArray(1:8).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedPos(13:18));
            [obj.dimerArray(1:6).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedPos(1:12));
            [obj.singleArray(1:12).Value] = deal(val{:});

            %load previous negative data
            obj.tripleArray(5).Value = CallingData.AddSelectedNeg(16);
            obj.doubleArray(9).Value = CallingData.AddSelectedNeg(15);
            val=num2cell(CallingData.AddSelectedNeg(11:14));
            [obj.dimerArray(7:10).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedNeg(1:10));
            [obj.singleArray(13:22).Value] = deal(val{:});
        end

        function obj = allSingleChanged(obj,~,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.singleArray,2));
            [obj.singleArray(:).Value] = val{:};
        end
        function obj = allDimerChanged(obj,~,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.dimerArray,2));
            [obj.dimerArray(:).Value] = val{:};
        end
        function obj = allTriplesChanged(obj,~,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.tripleArray,2));
            [obj.tripleArray(:).Value] = val{:};
        end
        function obj = allDualsChanged(obj,~,event)
            val = {event.Value};
            val = repmat(val,1,size(obj.doubleArray,2));
            [obj.doubleArray(:).Value] = val{:};
        end
    end
end