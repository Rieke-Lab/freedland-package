% Adapted from MHT's 'area summation figure'

classdef contrastReversingFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        preTime
        stimTime
        temporalFrequency
        monitorSampleRate
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
            ip.addParameter('monitorSampleRate', [], @(x)isvector(x));
            ip.parse(varargin{:});
            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            obj.temporalFrequency = ip.Results.temporalFrequency;
            obj.monitorSampleRate = ip.Results.monitorSampleRate;
            
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
            
            % Isolate spikes during stimTime
            epochResponseTrace = epochResponseTrace(prePts:(prePts+stimPts));
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            
            % Reinsert into binary vector
            epochResponseTrace = zeros(1,prePts:(prePts+stimPts));
            epochResponseTrace(S.sp) = 1;

            % Courtesy of M. Turner.
            L = length(epochResponseTrace); %length of signal, datapoints
            X = abs(fft(epochResponseTrace));
            X = X(1:L/2);
            f = obj.monitorSampleRate*(0:L/2-1)/L; %freq - hz
            [~, F1ind] = min(abs(f-obj.temporalFrequency)); %find index of F1 and F2 frequencies
            [~, F2ind] = min(abs(f-2*obj.temporalFrequency));
            F1 = 2*X(F1ind); %pA^2/Hz for current rec, (spikes/sec)^2/Hz for spike rate
            F2 = 2*X(F2ind); %double b/c of symmetry about zero
            
            % Pull associated parameter
            currentSpotSize = epoch.parameters('slices');
            obj.allSpotSizes = cat(1,obj.allSpotSizes,currentSpotSize);
            obj.allFirstPeriod = cat(1,obj.allFirstPeriod,F1);
            obj.allSecondPeriod = cat(1,obj.allSecondPeriod,F2);
            
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

