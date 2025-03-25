classdef NegAdducts < AddNLoss
%% Object for storing Negative Adduct Filter UI elements
    properties
        singleArray matlab.ui.control.CheckBox
        allSingles matlab.ui.control.CheckBox
        doubleArray matlab.ui.control.CheckBox
        dimerArray matlab.ui.control.CheckBox
        allDimers matlab.ui.control.CheckBox
        tripleArray matlab.ui.control.CheckBox
    end

    methods
        function obj = NegAdducts(CallingData)
            obj@AddNLoss(CallingData)

            %Construct negative adduct specific tick boxes
            %% Single Charged Adducts
            % Create S1
            obj.singleArray(1) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(1).Text = 'M+Cl';
            obj.singleArray(1).Position = [15 230 50 22];

            % Create S2
            obj.singleArray(2) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(2).Text = 'M+Br';
            obj.singleArray(2).Position = [15 210 51 22];

            % Create S3
            obj.singleArray(3) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(3).Text = 'M+Fac-H';
            obj.singleArray(3).Position = [15 190 72 22];

            % Create S4
            obj.singleArray(4) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(4).Text = 'M+HAc-H';
            obj.singleArray(4).Position = [15 170 74 22];

            % Create S5
            obj.singleArray(5) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(5).Text = 'M+NaAc-H';
            obj.singleArray(5).Position = [15 150 81 22];

            % Create S6
            obj.singleArray(6) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(6).Text = 'M+NH4Ac-H';
            obj.singleArray(6).Position = [15 130 81 22];

            % Create S7
            obj.singleArray(7) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(7).Text = 'M+TFA-H';
            obj.singleArray(7).Position = [15 110 74 22];

            % Create S8
            obj.singleArray(8) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(8).Text = 'M+Na-2H';
            obj.singleArray(8).Position = [15 90 74 22];

            % Create S9
            obj.singleArray(9) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(9).Text = 'M+Na-2H';
            obj.singleArray(9).Position = [15 70 74 22];

            % Create S10
            obj.singleArray(10) = uicheckbox(obj.SingleSubPanel);
            obj.singleArray(10).Text = 'M+K-2H';
            obj.singleArray(10).Position = [15 50 66 22];

            % Create AllSingle
            obj.allSingles = uicheckbox(obj.SingleSubPanel);
            obj.allSingles.ValueChangedFcn = @(src,event) {allSingleChanged(obj,src,event)};
            obj.allSingles.Text = 'All';
            obj.allSingles.Position = [15 250 35 22];

            %% Dimers
            % Create D1
            obj.dimerArray(1) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(1).Text = '2M-H';
            obj.dimerArray(1).Position = [15 90 51 22];

            % Create D2
            obj.dimerArray(2) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(2).Text = '2M+Fac-H';
            obj.dimerArray(2).Position = [15 70 78 22];

            % Create D3
            obj.dimerArray(3) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(3).Text = '2M+HAc-H';
            obj.dimerArray(3).Position = [15 50 81 22];

            % Create D4
            obj.dimerArray(4) = uicheckbox(obj.DimerSubPanel);
            obj.dimerArray(4).Text = '3M-H';
            obj.dimerArray(4).Position = [16 29 51 22];

            % Create AllD
            obj.allDimers = uicheckbox(obj.DimerSubPanel);
            obj.allDimers.ValueChangedFcn = @(src,event) {allDimerChanged(obj,src,event)};
            obj.allDimers.Text = 'All';
            obj.allDimers.Position = [15 110 35 22];

            %% Triple Charged
            % Create T1
            obj.tripleArray = uicheckbox(obj.TripleSubPanel);
            obj.tripleArray.Text = 'M-3H';
            obj.tripleArray.Position = [15 110 51 22];
            %% Dual Charged Adducts
            % Create Dual1
            obj.doubleArray = uicheckbox(obj.DualSubPanel);
            obj.doubleArray.Text = 'M-2H';
            obj.doubleArray.Position = [17 250 51 22];

            %load previous data
            obj.tripleArray(:).Value = CallingData.AddSelectedNeg(16);
            obj.doubleArray(:).Value = CallingData.AddSelectedNeg(15);
            val = num2cell(CallingData.AddSelectedNeg(11:14));
            [obj.dimerArray(:).Value] = deal(val{:});
            val = num2cell(CallingData.AddSelectedNeg(1:10));
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
    end
end