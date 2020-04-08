% Adapted from MHT's 'area summation figure'

classdef contrastDiskSizingFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        preTime
        stimTime
        type
    end
    
    properties (Access = private)
        axesHandle
        lineHandle
        fitLineHandle
        allFirstPeriod
        allSecondPeriod
        baselines
        allSpotSizes
        summaryData
    end
    
    methods
        
        function obj = contrastDiskSizingFigure(ampDevice, varargin)
            obj.ampDevice = ampDevice;            
            ip = inputParser();
            ip.addParameter('preTime', [], @(x)isvector(x));
            ip.addParameter('stimTime', [], @(x)isvector(x));
            ip.addParameter('temporalFrequency', [], @(x)isvector(x));
            ip.parse(varargin{:});
            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            obj.temporalFrequency = ip.Results.temporalFrequency;
            
            obj.createUi();
        end
        
        function createUi(obj)
            import appbox.*;
            
            obj.axesHandle = axes( ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle, obj.type);
            ylabel(obj.axesHandle, 'total spikes');
            title(obj.axesHandle,'find intersection');
        end

        
        function handleEpoch(obj, epoch)
            
            % Load amp data
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            
            prePts = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            iterations = obj.stimTime/1000 .* obj.temporalFrequency;
            iterationPts = stimPts / iterations; % Each cycle
            
            epochResponseTrace = epochResponseTrace(prePts:prePts+stimPts);

            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            
            newEpochResponse = zeros(iterations,2);
            for a = 1:iterations
                range = (a - 1) * iterationPts + 1 : a*iterationPts;
                cycle = S.sp(range);
                
                halfway = round(iterationPts / 2);
                newEpochResponse(a,1) = sum(cycle(1:halfway)); % First half-period
                newEpochResponse(a,2) = sum(cycle(halfway:end)); % Second half-period
            end
            newEpochResponse = nanmean(newEpochResponse,1); % Average over cycle
            
            % Pull associated parameter
            radiiUM = epoch.parameters('radii_um');
            radiiDim = epoch.parameters('radii_dimension');
            currentSpotSize = radiiUM(radiiDim);
            
            obj.allSpotSizes = cat(1,obj.allSpotSizes,currentSpotSize);
            obj.allFirstPeriod = cat(1,obj.allFirstPeriod,newEpochResponse(1));
            obj.allSecondPeriod = cat(1,obj.allSecondPeriod,newEpochResponse(2));
            
            obj.summaryData.spotSizes = unique(obj.allSpotSizes);
            obj.summaryData.meanResponses = zeros(size(obj.summaryData.spotSizes));
            for SpotSizeIndex = 1:length(obj.summaryData.spotSizes)
                pullIndices = (obj.summaryData.spotSizes(SpotSizeIndex) == obj.allSpotSizes);
                obj.summaryData.meanResponsesFirst(SpotSizeIndex) = mean(obj.allFirstPeriod(pullIndices));
                obj.summaryData.meanResponsesSecond(SpotSizeIndex) = mean(obj.allSecondPeriod(pullIndices));
            end
            
            if isempty(obj.lineHandle)
                obj.lineHandle = line(obj.summaryData.spotSizes, obj.summaryData.meanResponsesFirst,...
                    'Parent', obj.axesHandle,'Color','r','Marker','o');
                obj.lineHandle = line(obj.summaryData.spotSizes, obj.summaryData.meanResponsesSecond,...
                    'Parent', obj.axesHandle,'Color','b','Marker','o');
            else
                set(obj.lineHandle, 'XData', obj.summaryData.spotSizes,...
                    'YData', obj.summaryData.meanResponsesFirst);
                set(obj.lineHandle, 'XData', obj.summaryData.spotSizes,...
                    'YData', obj.summaryData.meanResponsesSecond);
            end
        end
        
    end
    
end

