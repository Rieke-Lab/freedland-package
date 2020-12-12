% Adapted from MHT's 'area summation figure'

classdef contrastReversingFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        preTime
        stimTime
        temporalFrequency
    end
    
    properties (Access = private)
        axesHandle
        lineHandle
        lineHandle2
        fitLineHandle
        allFirstPeriod
        allSecondPeriod
        baselines
        allSpotSizes
        summaryData
    end
    
    methods
        
        function obj = contrastReversingFigure(ampDevice, varargin)
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
            iconDir = [fileparts(fileparts(mfilename('fullpath'))), '\+utils\+icons\'];
            toolbar = findall(obj.figureHandle, 'Type', 'uitoolbar');
            
            obj.axesHandle = axes( ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle, 'divisions');
            ylabel(obj.axesHandle, 'total spikes');
            
            intersectionButton = uipushtool( ...
                'Parent', toolbar, ...
                'TooltipString', 'intersection', ...
                'Separator', 'on', ...
                'ClickedCallback', @obj.intersection);
            setIconImage(intersectionButton, [iconDir, 'DoG.png']);
        end

        function handleEpoch(obj, epoch)
            
            % Load amp data
            response            = epoch.getResponse(obj.ampDevice);
            epochResponseTrace  = response.getData();
            sampleRate          = response.sampleRate.quantityInBaseUnits;
            
            % Define experiment parameters
            prePts  = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            iterations      = obj.stimTime/1000 .* obj.temporalFrequency;
            iterationPts    = stimPts / iterations; % Each cycle
            halfway = round(iterationPts / 2);
            
            % Isolate spikes during stimTime
            epochResponseTrace = epochResponseTrace(prePts:(prePts+stimPts));
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);

            % Identify phase-specific spikes
            newEpochResponse = zeros(iterations,2);
            for a = 1:iterations
                firstPt = ((a - 1) * iterationPts) + 1;
                newEpochResponse(a,1) = sum(S.sp > firstPt & S.sp < (firstPt + halfway));
                newEpochResponse(a,2) = sum(S.sp > (firstPt + halfway) & S.sp < (firstPt + iterationPts));
            end

            newEpochResponse = nanmean(newEpochResponse,1); % Average over cycle
            
            % Pull associated parameter
            currentSpotSize = epoch.parameters('slices');
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
                obj.lineHandle2 = line(obj.summaryData.spotSizes, obj.summaryData.meanResponsesSecond,...
                    'Parent', obj.axesHandle,'Color','b','Marker','o');
            else
                set(obj.lineHandle, 'XData', obj.summaryData.spotSizes,...
                    'YData', obj.summaryData.meanResponsesFirst);
                set(obj.lineHandle2, 'XData', obj.summaryData.spotSizes,...
                    'YData', obj.summaryData.meanResponsesSecond);
            end
        end
        
    end
    
     methods (Access = private)
        
        function intersection(obj, ~, ~)
            x = obj.lineHandle.XData;
            y1 = obj.lineHandle.YData;
            y2 = obj.lineHandle2.YData;
            
            resolution = 0.1; % um
            xq = min(x):resolution:max(x);
            
            y1_i = interp1(x,y1,xq);
            y2_i = interp1(x,y2,xq);

            crossZero = abs(diff(sign(y1_i - y2_i))) > 0;
            crossZero = unique(round(xq(crossZero)));
            
            title(obj.axesHandle,num2str(crossZero));
        end

    end
    
end

