classdef PosAdducts < AddNLoss
    %% Object for storing Positive adduct filter UI elements
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
        function obj = PosAdducts(CallingData)
            obj@AddNLoss(CallingData)

            %Construct Positive Mode specific Tick boxes
            %% Single Charged Adducts
            % Create S1
            obj.singleArray(1) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(1).Text = 'M+NH4';
            obj.singleArray(1).Position = [15 230 63 22];

            % Create S2
            obj.singleArray(2) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(2).Text = 'M+Na';
            obj.singleArray(2).Position = [15 210 54 22];

            % Create S3
            obj.singleArray(3) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(3).Text = 'M+K';
            obj.singleArray(3).Position = [15 190 47 22];

            % Create S4
            obj.singleArray(4) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(4).Text = 'M+ACN+H';
            obj.singleArray(4).Position = [15 170 80 22];

            % Create S5
            obj.singleArray(5) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(5).Text = 'M+MeOH+H';
            obj.singleArray(5).Position = [15 150 89 22];

            % Create S6
            obj.singleArray(6) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(6).Text = 'M+IPA+H';
            obj.singleArray(6).Position = [15 130 73 22];

            % Create S7
            obj.singleArray(7) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(7).Text = 'M+DMSO+H';
            obj.singleArray(7).Position = [15 110 91 22];

            % Create S8
            obj.singleArray(8) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(8).Text = 'M+ACN+Na';
            obj.singleArray(8).Position = [15 90 87 22];

            % Create S9
            obj.singleArray(9) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(9).Text = 'M+2ACN+H';
            obj.singleArray(9).Position = [15 70 87 22];

            % Create S10
            obj.singleArray(10) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(10).Text = 'M+IPA+Na';
            obj.singleArray(10).Position = [15 50 80 22];

            % Create S11
            obj.singleArray(11) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(11).Text = 'M+2Na-H';
            obj.singleArray(11).Position = [15 30 74 22];

            % Create S12
            obj.singleArray(12) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(12).Text = 'M+2K-H';
            obj.singleArray(12).Position = [15 10 66 22];

            % Create AllSingle
            obj.allSingles = uicheckbox(obj.SingleSubPanel);
            obj.allSingles.ValueChangedFcn = @(src,event) {allSingleChanged(obj,src,event)};
            obj.allSingles.Text = 'All';
            obj.allSingles.Position = [15 250 35 22];

            %% Dimers
            % Create D1
            obj.dimerArray(1) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(1).Text = '2M+H';
            obj.dimerArray(1).Position = [15 105 54 22];

            % Create D2
            obj.dimerArray(2) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(2).Text = '2M+NH4';
            obj.dimerArray(2).Position = [15 85 70 22];

            % Create D3
            obj.dimerArray(3) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(3).Text = '2M+Na';
            obj.dimerArray(3).Position = [15 65 61 22];

            % Create D4
            obj.dimerArray(4) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(4).Text = '2M+K';
            obj.dimerArray(4).Position = [15 45 54 22];

            % Create D5
            obj.dimerArray(5) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(5).Text = '2M+ACN+H';
            obj.dimerArray(5).Position = [15 25 87 22];

            % Create D6
            obj.dimerArray(6) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(6).Text = '2M+ACN+Na';
            obj.dimerArray(6).Position = [15 5 93 22];

            % Create AllD
            obj.allDimers = uicheckbox(obj.DimerSubPanel);
            obj.allDimers.ValueChangedFcn = @(src,event) {allDimerChanged(obj,src,event)};
            obj.allDimers.Text = 'All';
            obj.allDimers.Position = [15 125 35 22];

            %% Triple Charged
            % Create T1
            obj.tripleArray(1) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(1).Text = 'M+3H';
            obj.tripleArray(1).Position = [15 105 54 22];

            % Create T2
            obj.tripleArray(2) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(2).Text = 'M+Na+2H';
            obj.tripleArray(2).Position = [15 85 77 22];

            % Create T3
            obj.tripleArray(3) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(3).Text = 'M+2Na+H';
            obj.tripleArray(3).Position = [15 65 77 22];

            % Create T4
            obj.tripleArray(4) = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray(4).Text = 'M+3Na';
            obj.tripleArray(4).Position = [15 45 61 22];

            % Create All Triples
            obj.allTriples = uicheckbox(obj.TripleSubPanel);
            obj.allTriples.ValueChangedFcn = @(src,event) {allTriplesChanged(obj,src,event)};
            obj.allTriples.Text = 'All';
            obj.allTriples.Position = [15 125 35 22];
            %% Dual Charged Adducts
            % Create Dual1
            obj.doubleArray(1) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(1).Text = 'M+2H';
            obj.doubleArray(1).Position = [15 230 54 22];

            % Create Dual2
            obj.doubleArray(2) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(2).Text = 'M+NH4+H';
            obj.doubleArray(2).Position = [15 210 79 22];

            % Create Dual3
            obj.doubleArray(3) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(3).Text = 'M+Na+H';
            obj.doubleArray(3).Position = [15 190 70 22];

            % Create Dual4
            obj.doubleArray(4) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(4).Text = 'M+K+H';
            obj.doubleArray(4).Position = [15 170 63 22];

            % Create Dual5
            obj.doubleArray(5) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(5).Text = 'M+2Na';
            obj.doubleArray(5).Position = [15 150 61 22];

            % Create Dual6
            obj.doubleArray(6) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(6).Text = 'M+ACN+2H';
            obj.doubleArray(6).Position = [15 130 87 22];

            % Create Dual7
            obj.doubleArray(7) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(7).Text = 'M+2ACN+2H';
            obj.doubleArray(7).Position = [15 110 93 22];

            % Create Dual8
            obj.doubleArray(8) = uicheckbox(obj.DualSubPanel);
            obj.doubleArray(8).Text = 'M+3ACN+2H';
            obj.doubleArray(8).Position = [15 90 93 22];
            % Create All Dual
            obj.allDoubles = uicheckbox(obj.DualSubPanel);
            obj.allDoubles.ValueChangedFcn = @(src,event) {allDualsChanged(obj,src,event)};
            obj.allDoubles.Text = 'All';
            obj.allDoubles.Position = [15 250 35 22];

            %load previous data
            val = num2cell(CallingData.AddSelectedPos(27:30));
            [obj.tripleArray(:).Value] = deal(val{:});
            val = num2cell(CallingData.AddSelectedPos(19:26));
            [obj.doubleArray(:).Value] = deal(val{:});
            val = num2cell(CallingData.AddSelectedPos(13:18));
            [obj.dimerArray(:).Value] = deal(val{:});
            val = num2cell(CallingData.AddSelectedPos(1:12));
            [obj.singleArray(:).Value] = deal(val{:});
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