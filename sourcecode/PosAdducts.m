classdef PosAdducts < AddNLoss
    % BSD 3-Clause License
% 
% Copyright (c) 2022, Adrian Haun
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
        function obj = PosAdducts(CallingData)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj@AddNLoss(CallingData)

            %Construct Pos Specific Tickboxes
            %% Single Charged Adducts
            % Create S1
            obj.SingleArray(1) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(1).Text = 'M+NH4';
            obj.SingleArray(1).Position = [15 230 63 22];

            % Create S2
            obj.SingleArray(2) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(2).Text = 'M+Na';
            obj.SingleArray(2).Position = [15 210 54 22];

            % Create S3
            obj.SingleArray(3) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(3).Text = 'M+K';
            obj.SingleArray(3).Position = [15 190 47 22];

            % Create S4
            obj.SingleArray(4) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(4).Text = 'M+ACN+H';
            obj.SingleArray(4).Position = [15 170 80 22];

            % Create S5
            obj.SingleArray(5) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(5).Text = 'M+MeOH+H';
            obj.SingleArray(5).Position = [15 150 89 22];

            % Create S6
            obj.SingleArray(6) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(6).Text = 'M+IPA+H';
            obj.SingleArray(6).Position = [15 130 73 22];

            % Create S7
            obj.SingleArray(7) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(7).Text = 'M+DMSO+H';
            obj.SingleArray(7).Position = [15 110 91 22];

            % Create S8
            obj.SingleArray(8) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(8).Text = 'M+ACN+Na';
            obj.SingleArray(8).Position = [15 90 87 22];

            % Create S9
            obj.SingleArray(9) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(9).Text = 'M+2ACN+H';
            obj.SingleArray(9).Position = [15 70 87 22];

            % Create S10
            obj.SingleArray(10) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(10).Text = 'M+IPA+Na';
            obj.SingleArray(10).Position = [15 50 80 22];

            % Create S11
            obj.SingleArray(11) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(11).Text = 'M+2Na-H';
            obj.SingleArray(11).Position = [15 30 74 22];

            % Create S12
            obj.SingleArray(12) = uicheckbox(obj.SingleSubPanel);
            obj.SingleArray(12).Text = 'M+2K-H';
            obj.SingleArray(12).Position = [15 10 66 22];

            % Create AllSingle
            obj.AllSingles = uicheckbox(obj.SingleSubPanel);
            obj.AllSingles.ValueChangedFcn = @(src,event) {AllSingleChanged(obj,src,event)};
            obj.AllSingles.Text = 'All';
            obj.AllSingles.Position = [15 250 35 22];

            %% Dimers
            % Create D1
            obj.DimerArray(1) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(1).Text = '2M+H';
            obj.DimerArray(1).Position = [15 105 54 22];

            % Create D2
            obj.DimerArray(2) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(2).Text = '2M+NH4';
            obj.DimerArray(2).Position = [15 85 70 22];

            % Create D3
            obj.DimerArray(3) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(3).Text = '2M+Na';
            obj.DimerArray(3).Position = [15 65 61 22];

            % Create D4
            obj.DimerArray(4) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(4).Text = '2M+K';
            obj.DimerArray(4).Position = [15 45 54 22];

            % Create D5
            obj.DimerArray(5) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(5).Text = '2M+ACN+H';
            obj.DimerArray(5).Position = [15 25 87 22];

            % Create D6
            obj.DimerArray(6) = uicheckbox(obj.DimerSubPanel);
            obj.DimerArray(6).Text = '2M+ACN+Na';
            obj.DimerArray(6).Position = [15 5 93 22];

            % Create AllD
            obj.AllDimers = uicheckbox(obj.DimerSubPanel);
            obj.AllDimers.ValueChangedFcn = @(src,event) {AllDimerChanged(obj,src,event)};
            obj.AllDimers.Text = 'All';
            obj.AllDimers.Position = [15 125 35 22];

            %% Triple Charged
            % Create T1
            obj.TripleArray(1) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(1).Text = 'M+3H';
            obj.TripleArray(1).Position = [15 105 54 22];

            % Create T2
            obj.TripleArray(2) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(2).Text = 'M+Na+2H';
            obj.TripleArray(2).Position = [15 85 77 22];

            % Create T3
            obj.TripleArray(3) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(3).Text = 'M+2Na+H';
            obj.TripleArray(3).Position = [15 65 77 22];

            % Create T4
            obj.TripleArray(4) = uicheckbox(obj.TripleSubPanel);
            obj.TripleArray(4).Text = 'M+3Na';
            obj.TripleArray(4).Position = [15 45 61 22];

            % Create All Triples
            obj.AllTriples = uicheckbox(obj.TripleSubPanel);
            obj.AllTriples.ValueChangedFcn = @(src,event) {AllTriplesChanged(obj,src,event)};
            obj.AllTriples.Text = 'All';
            obj.AllTriples.Position = [15 125 35 22];
            %% Dual Charged Adducts
            % Create Dual1
            obj.DoubleArray(1) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(1).Text = 'M+2H';
            obj.DoubleArray(1).Position = [15 230 54 22];

            % Create Dual2
            obj.DoubleArray(2) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(2).Text = 'M+NH4+H';
            obj.DoubleArray(2).Position = [15 210 79 22];

            % Create Dual3
            obj.DoubleArray(3) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(3).Text = 'M+Na+H';
            obj.DoubleArray(3).Position = [15 190 70 22];

            % Create Dual4
            obj.DoubleArray(4) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(4).Text = 'M+K+H';
            obj.DoubleArray(4).Position = [15 170 63 22];

            % Create Dual5
            obj.DoubleArray(5) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(5).Text = 'M+2Na';
            obj.DoubleArray(5).Position = [15 150 61 22];

            % Create Dual6
            obj.DoubleArray(6) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(6).Text = 'M+ACN+2H';
            obj.DoubleArray(6).Position = [15 130 87 22];

            % Create Dual7
            obj.DoubleArray(7) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(7).Text = 'M+2ACN+2H';
            obj.DoubleArray(7).Position = [15 110 93 22];

            % Create Dual8
            obj.DoubleArray(8) = uicheckbox(obj.DualSubPanel);
            obj.DoubleArray(8).Text = 'M+3ACN+2H';
            obj.DoubleArray(8).Position = [15 90 93 22];
            % Create All Dual
            obj.AllDoubles = uicheckbox(obj.DualSubPanel);
            obj.AllDoubles.ValueChangedFcn = @(src,event) {AllDualsChanged(obj,src,event)};
            obj.AllDoubles.Text = 'All';
            obj.AllDoubles.Position = [15 250 35 22];
            %load previous data
            val = num2cell(CallingData.AddSelectedPos(27:30));
            [obj.TripleArray(:).Value] = deal(val{:});
            val = num2cell(CallingData.AddSelectedPos(19:26));
            [obj.DoubleArray(:).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedPos(13:18));
            [obj.DimerArray(:).Value] = deal(val{:});
            val=num2cell(CallingData.AddSelectedPos(1:12));
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
        function obj = AllTriplesChanged(obj,src,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.TripleArray,2));
            [obj.TripleArray(:).Value] = val{:};
        end
        function obj = AllDualsChanged(obj,src,event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = {event.Value};
            val = repmat(val,1,size(obj.DoubleArray,2));
            [obj.DoubleArray(:).Value] = val{:};
        end
    end
end