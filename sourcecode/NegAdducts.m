classdef NegAdducts < AddNLoss
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        SingleArray matlab.ui.control.CheckBox
        AllSingles matlab.ui.control.CheckBox
        DoubleArray matlab.ui.control.CheckBox
        DimerArray matlab.ui.control.CheckBox
        AllDimers matlab.ui.control.CheckBox
        TripleArray matlab.ui.control.CheckBox
    end

    methods
        function obj = NegAdducts(CallingData)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj@AddNLoss(CallingData)

            %Construct Neg Specific Tickboxes
            %% Single Charged Adducts
            % Create S1
            obj.SingleArray(1) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(1).Text = 'M+Cl';
            obj.SingleArray(1).Position = [15 230 50 22];

            % Create S2
            obj.SingleArray(2) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(2).Text = 'M+Br';
            obj.SingleArray(2).Position = [15 210 51 22];

            % Create S3
            obj.SingleArray(3) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(3).Text = 'M+Fac-H';
            obj.SingleArray(3).Position = [15 190 72 22];

            % Create S4
            obj.SingleArray(4) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(4).Text = 'M+HAc-H';
            obj.SingleArray(4).Position = [15 170 74 22];

            % Create S5
            obj.SingleArray(5) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(5).Text = 'M+NaAc-H';
            obj.SingleArray(5).Position = [15 150 81 22];

            % Create S6
            obj.SingleArray(6) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(6).Text = 'M+NH4Ac-H';
            obj.SingleArray(6).Position = [15 130 81 22];

            % Create S7
            obj.SingleArray(7) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(7).Text = 'M+TFA-H';
            obj.SingleArray(7).Position = [15 110 74 22];

            % Create S8
            obj.SingleArray(8) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(8).Text = 'M+Na-2H';
            obj.SingleArray(8).Position = [15 90 74 22];

            % Create S9
            obj.SingleArray(9) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(9).Text = 'M+Na-2H';
            obj.SingleArray(9).Position = [15 70 74 22];

            % Create S10
            obj.SingleArray(10) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(10).Text = 'M+K-2H';
            obj.SingleArray(10).Position = [15 50 66 22];

            % Create AllSingle
            obj.AllSingles = uicheckbox(obj.SingleSubPanel);
            obj.AllSingles.ValueChangedFcn = @(src,event) {AllSingleChanged(obj,src,event)};
            obj.AllSingles.Text = 'All';
            obj.AllSingles.Position = [15 250 35 22];

            %% Dimers
            % Create D1
            obj.DimerArray(1) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(1).Text = '2M-H';
            obj.DimerArray(1).Position = [15 90 51 22];

            % Create D2
            obj.DimerArray(2) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(2).Text = '2M+Fac-H';
            obj.DimerArray(2).Position = [15 70 78 22];

            % Create D3
            obj.DimerArray(3) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(3).Text = '2M+HAc-H';
            obj.DimerArray(3).Position = [15 50 81 22];

            % Create D4
            obj.DimerArray(4) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(4).Text = '3M-H';
            obj.DimerArray(4).Position = [16 29 51 22];

            % Create AllD
            obj.AllDimers = uicheckbox(obj.DimerSubPanel);
            obj.AllDimers.ValueChangedFcn = @(src,event) {AllDimerChanged(obj,src,event)};
            obj.AllDimers.Text = 'All';
            obj.AllDimers.Position = [15 110 35 22];

            %% Triple Charged
            % Create T1
            obj.TripleArray = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray.Text = 'M-3H';
            obj.TripleArray.Position = [15 110 51 22];
            %% Dual Charged Adducts
            % Create Dual1
            obj.DoubleArray = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray.Text = 'M-2H';
            obj.DoubleArray.Position = [17 250 51 22];

            %load previous data
            obj.TripleArray(:).Value = CallingData.AddSelectedNeg(16);
            obj.DoubleArray(:).Value = CallingData.AddSelectedNeg(15);
            val=num2cell(CallingData.AddSelectedNeg(11:14));
            [obj.DimerArray(:).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedNeg(1:10));
            [obj.SingleArray(:).Value] = deal(val{:});
        end

        function obj = AllSingleChanged(obj,src,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.SingleArray,2));
            [obj.SingleArray(:).Value] = val{:};
        end
        function obj = AllDimerChanged(obj,src,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.DimerArray,2));
            [obj.DimerArray(:).Value] = val{:};
        end
    end
end