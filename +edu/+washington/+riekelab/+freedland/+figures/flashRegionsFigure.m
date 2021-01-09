% Adapted from MHT's 'area summation figure'

classdef flashRegionsFigure < symphonyui.core.FigureHandler
    
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
        allEpochResponses
        baselines
        allSpotSizes
        summaryData
        rawExport
        allTrackingSpotSize
    end
    
    methods
        
        function obj = flashRegionsFigure(ampDevice, varargin)
            obj.ampDevice = ampDevice;            
            ip = inputParser();
            ip.addParameter('preTime', [], @(x)isvector(x));
            ip.addParameter('stimTime', [], @(x)isvector(x));
            ip.addParameter('type', []);
            ip.parse(varargin{:});
            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            obj.type = ip.Results.type;
            
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
            title(obj.axesHandle,'rf filter fit (onset - offset)');
        end

        
        function handleEpoch(obj, epoch)
            
            % Load amp data
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            currentSpotSize = sum(epoch.parameters('flashedRegions'));
            trackingSpotSize = epoch.parameters('flashedRegions');
            
            prePts = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            preScaleFactor = stimPts / prePts;
            
            % Count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            newEpochResponse = sum(S.sp > prePts & S.sp <= prePts + stimPts) -... % spike count during onset
                sum(S.sp >= prePts + stimPts);                                    % spike count during offset
            newBaseline = preScaleFactor * sum(S.sp <= prePts); % spike count before stim, scaled by length
         
            obj.allSpotSizes = cat(1,obj.allSpotSizes,currentSpotSize);
            obj.allTrackingSpotSize = cat(1,obj.allTrackingSpotSize,trackingSpotSize);
            obj.allEpochResponses = cat(1,obj.allEpochResponses,newEpochResponse);
            obj.baselines = cat(1,obj.baselines,newBaseline);
            
            % For graphing
            obj.summaryData.spotSizes = unique(obj.allSpotSizes);
            obj.summaryData.meanResponses = zeros(size(obj.summaryData.spotSizes));
            for SpotSizeIndex = 1:length(obj.summaryData.spotSizes)
                pullIndices = (obj.summaryData.spotSizes(SpotSizeIndex) == obj.allSpotSizes);
                obj.summaryData.meanResponses(SpotSizeIndex) = mean(obj.allEpochResponses(pullIndices));
            end
            
            % For calculations
            obj.summaryData.trackingSpotSizes = unique(obj.allTrackingSpotSize,'rows');
            obj.summaryData.trackingResponses = zeros(size(obj.summaryData.trackingSpotSizes,1),1);
            for SpotSizeIndex = 1:size(obj.summaryData.trackingSpotSizes,1)
                [~,pullIndices] = ismember(obj.summaryData.trackingSpotSizes(SpotSizeIndex,:),obj.allTrackingSpotSize,'rows');
                obj.summaryData.trackingResponses(SpotSizeIndex) = mean(obj.allEpochResponses(pullIndices));
            end
            
            % Plot
            if isempty(obj.lineHandle)
                obj.lineHandle = line(obj.summaryData.spotSizes, obj.summaryData.meanResponses,...
                    'Parent', obj.axesHandle,'Color','k','Marker','o');
            else
                set(obj.lineHandle, 'XData', obj.summaryData.spotSizes,...
                    'YData', obj.summaryData.meanResponses);
            end
            
            %%%%% Save data for further experiments
            if rank(obj.summaryData.trackingSpotSizes) == size(obj.summaryData.trackingSpotSizes,2)
                
                % A * w = B, where w is weights
                w = obj.summaryData.trackingSpotSizes\obj.summaryData.trackingResponses;

                % Export
                title(obj.axesHandle,num2str(w));
                w = [w, w/max(w)]'; % Add row for normalization
                dlmwrite(strcat('Documents/weights.txt'),w)
            end
            %%%%%
        end
    end
    
end

