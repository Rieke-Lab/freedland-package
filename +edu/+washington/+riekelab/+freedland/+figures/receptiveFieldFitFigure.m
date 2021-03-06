% Adapted from MHT's 'area summation figure'

classdef receptiveFieldFitFigure < symphonyui.core.FigureHandler
    
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
        allCenterR
        allAnnulusR
        allSurroundR
    end
    
    methods
        
        function obj = receptiveFieldFitFigure(ampDevice, varargin)
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
            
            %load amp data
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;

            if strfind(obj.type,'center') > 0
                currentSpotSize = epoch.parameters('rfSigmaCenter');
            elseif strfind(obj.type,'surround')
                currentSpotSize = epoch.parameters('rfSigmaSurround');
            elseif strfind(obj.type,'rotation')
                currentSpotSize = epoch.parameters('rotation');
            elseif strfind(obj.type,'experimentID')
                currentSpotSize = epoch.parameters('experimentID');
            end
            
            prePts = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            preScaleFactor = stimPts / prePts;
            
            %count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            newEpochResponse = sum(S.sp > prePts & S.sp <= prePts + stimPts) -... % spike count during onset
                sum(S.sp >= prePts + stimPts);                                    % spike count during offset
            newBaseline = preScaleFactor * sum(S.sp <= prePts); % spike count before stim, scaled by length
         
            obj.allSpotSizes = cat(1,obj.allSpotSizes,currentSpotSize);
            obj.allEpochResponses = cat(1,obj.allEpochResponses,newEpochResponse);
            obj.baselines = cat(1,obj.baselines,newBaseline);
            
            obj.summaryData.spotSizes = unique(obj.allSpotSizes);
            obj.summaryData.meanResponses = zeros(size(obj.summaryData.spotSizes));
            for SpotSizeIndex = 1:length(obj.summaryData.spotSizes)
                pullIndices = (obj.summaryData.spotSizes(SpotSizeIndex) == obj.allSpotSizes);
                obj.summaryData.meanResponses(SpotSizeIndex) = mean(obj.allEpochResponses(pullIndices));
            end
            
            if isempty(obj.lineHandle)
                obj.lineHandle = line(obj.summaryData.spotSizes, obj.summaryData.meanResponses,...
                    'Parent', obj.axesHandle,'Color','k','Marker','o');
            else
                set(obj.lineHandle, 'XData', obj.summaryData.spotSizes,...
                    'YData', obj.summaryData.meanResponses);
            end
        end
        
    end
    
end

