classdef GCData < RawData
    % Class for storing group settings and performing functions from Raw
    % data until Feature data stage
    properties
        separationType (1,1) string = "GC"
    end

    methods
        function obj = GCData(groupName,window)
            %Construct an instance of this class
            if nargin == 0
                groupName = 0;
                window = 0;
            end
            obj = obj@RawData(groupName,window);
            % set default parameters
            obj.withinFileMassTolerance = 0.1;
            obj.withinFileMassUnit = "Da";
            obj.roiMinOccurrence = 10;
            obj.peakMinWidth = 0.8;
            obj.peakMaxWidth = 10;
            obj.betweenFileMassTolerance = 0.05;
            obj.betweenFileMassUnit = "Da";
        end

        function [Output,obj] = processGroup(obj,varargin)
            [Output,obj] = obj.extractFeaturesFromMassData(varargin);
        end

    end
end