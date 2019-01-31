% Courtesy of Max Turner
classdef AreaSummationNormFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        recordingType
        preTime
        stimTime
    end
    
    properties (Access = private)
        axesHandle
        lineHandle
        fitLineHandle
        allEpochResponses
        baselines
        allSpotSizes
        summaryData
    end
    
    methods
        
        function obj = AreaSummationFigure(ampDevice, varargin)
            obj.ampDevice = ampDevice;            
            ip = inputParser();
            ip.addParameter('recordingType', [], @(x)ischar(x));
            ip.addParameter('preTime', [], @(x)isvector(x));
            ip.addParameter('stimTime', [], @(x)isvector(x));
            ip.parse(varargin{:});
            obj.recordingType = ip.Results.recordingType;
            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            
            obj.createUi();
        end
        
        function createUi(obj)

            import appbox.*;
            iconDir = '+edu/+washington/+riekelab/+freedland/+utils/+icons/';
            toolbar = findall(obj.figureHandle, 'Type', 'uitoolbar');
            
            fitDoGButton = uipushtool( ...
                'Parent', toolbar, ...
                'TooltipString', 'Fit DoG', ...
                'Separator', 'on', ...
                'ClickedCallback', @obj.onSelectedFitDoG);
            setIconImage(fitDoGButton, [iconDir, 'DoG.png']);
            
            
            obj.axesHandle = axes( ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle, 'Spot Diameter (um)');
            ylabel(obj.axesHandle, 'Response');
            title(obj.axesHandle,'Area summation curve');
            
        end

        
        function handleEpoch(obj, epoch)
            %load amp data
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            currentSpotSize = epoch.parameters('currentSpotSize');
            prePts = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            preScaleFactor = stimPts / prePts;
            
            epochResponseTrace = epochResponseTrace(1:prePts+stimPts);
            %count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            newEpochResponse = sum(S.sp > prePts); %spike count during stim
            newBaseline = preScaleFactor * sum(S.sp < prePts); %spike count before stim, scaled by length

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
    
    methods (Access = private)
        
        function onSelectedFitDoG(obj, ~, ~)
            % normalize -- by JMF
            obj.summaryData.meanResponses = obj.summaryData.meanResponses - obj.summaryData.meanResponses(1,1);
            
            params0 = [max(obj.summaryData.meanResponses) / 2,50,...
                max(obj.summaryData.meanResponses) / 2, 150];
            [Kc, sigmaC, Ks, sigmaS] = ...
                edu.washington.riekelab.freedland.utils.fitDoGAreaSummation(obj.summaryData.spotSizes,obj.summaryData.meanResponses,params0);
            fitX = 0:(1.1*max(obj.summaryData.spotSizes));
            fitY = edu.washington.riekelab.freedland.utils.DoGAreaSummation([Kc sigmaC Ks sigmaS], fitX);
            
            if isempty(obj.fitLineHandle)
                obj.fitLineHandle = line(fitX, fitY, 'Parent', obj.axesHandle);
            else
                set(obj.fitLineHandle, 'XData', fitX,...
                    'YData', fitY);
            end
            set(obj.fitLineHandle,'Color',[1 0 0],'LineWidth',2,'Marker','none');
            tempKc = Kc / (Kc + Ks);
            str = {['SigmaC = ',num2str(sigmaC)],['sigmaS = ',num2str(sigmaS)],...
            ['Kc = ',num2str(tempKc)]};
            title(obj.axesHandle,str);
        end

    end
    
end

