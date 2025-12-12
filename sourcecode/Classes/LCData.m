classdef LCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        separationType (1,1) string = "LC"
    end

    methods
        function obj = LCData(groupName,window)
            %Construct an instance of this class
            if nargin == 0
                groupName = 0;
                window = 0;
            end
            obj = obj@RawData(groupName,window);
        end

        function [Output,obj] = processGroup(obj,varargin)
            [Output,obj] = obj.extractFeaturesFromMassData(varargin);
        end
        
    end
end