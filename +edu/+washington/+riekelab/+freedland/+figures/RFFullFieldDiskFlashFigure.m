% Adapted from MHT's 'area summation figure'

classdef RFFullFieldDiskFlashFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        preTime
        stimTime
        tailTime
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
        
        function obj = RFFullFieldDiskFlashFigure(ampDevice, varargin)
            obj.ampDevice = ampDevice;            
            ip = inputParser();
            ip.addParameter('preTime', [], @(x)isvector(x));
            ip.addParameter('stimTime', [], @(x)isvector(x));
            ip.addParameter('tailTime', [], @(x)isvector(x));
            ip.parse(varargin{:});

            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            obj.tailTime = ip.Results.tailTime;
            
            obj.createUi();
        end
        
        function createUi(obj)
            import appbox.*;
            
            obj.axesHandle = axes( ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle, 'preset number');
            ylabel(obj.axesHandle, 'ON flashed disk spikes / ON image spikes');
            title(obj.axesHandle,'full field flash quality');
            
        end

        
        function handleEpoch(obj, epoch)
            
            %load amp data
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            currentSpotSize = epoch.parameters('presetNo');
            prePts = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            tailPts = sampleRate*obj.tailTime/1000;
            
            epochResponseTrace = epochResponseTrace(1:2*(prePts+stimPts+tailPts));
            totalStimTime = prePts+stimPts+tailPts;
            
            %count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            
            % pull total spikes during image flashes
            naturalImageRespON = sum(S.sp > prePts & S.sp < prePts + stimPts); %spike count during stim
            naturalImageRespOFF = sum(S.sp > prePts + stimPts & S.sp < prePts + stimPts + tailPts); %spike count during stim
            flashRespON = sum(S.sp > totalStimTime + prePts & S.sp < totalStimTime + prePts + stimPts); %spike count during stim
            flashRespOFF = sum(S.sp > totalStimTime + prePts + stimPts & S.sp < totalStimTime + prePts + stimPts + tailPts); %spike count during stim
            newEpochResponse = (flashRespON-flashRespOFF) / (naturalImageRespON-naturalImageRespOFF);
            
            if newEpochResponse == Inf || newEpochResponse == -Inf || isnan(newEpochResponse)
                newEpochResponse = 0;
            end
            
            obj.allSpotSizes = cat(1,obj.allSpotSizes,currentSpotSize);
            obj.allEpochResponses = cat(1,obj.allEpochResponses,newEpochResponse);
            
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

