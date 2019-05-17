% Makes three plots: a radial map, DoG RF, and an empirical
% receptive field.

classdef RFRadialMapFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        spinningFrequency
        expandingFrequency
        frameMonitor
        stageDevice
        canvasSize
        ovalWidth
        umPerPix
    end
    
    properties (Access = private)
        axesHandle
        lineHandle
        fitLineHandle
        radialData
        RFMapData
        spins
        growth
        summaryData
        allRadii
        allAngles
        allResponses
        allRFResponses
        totalMeshGrid
    end
    
    methods
        
        
        function obj = RFRadialMapFigure(ampDevice, stage, frameMonitor, varargin)
            obj.ampDevice = ampDevice;
            obj.stageDevice = stage;
            obj.frameMonitor = frameMonitor;
            ip = inputParser();
            ip.addParameter('spinningFrequency', 2);
            ip.addParameter('expandingFrequency', 100);
            ip.addParameter('canvasSize', [800 600]);
            ip.addParameter('ovalWidth', 5);
            ip.addParameter('umPerPix', 1.2);
            ip.parse(varargin{:});
            obj.spinningFrequency = ip.Results.spinningFrequency;
            obj.expandingFrequency = ip.Results.expandingFrequency;
            obj.canvasSize = ip.Results.canvasSize;
            obj.ovalWidth = ip.Results.ovalWidth;
            obj.umPerPix = ip.Results.umPerPix;
            obj.createUi();
        end
        
        function createUi(obj)
            import appbox.*;

            obj.axesHandle(1) = subplot(1,2,1, ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle(1), 'orientation');
            ylabel(obj.axesHandle(1), 'response');
            
            obj.axesHandle(2) = subplot(1,2,2, ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle(2), 'distance from center (um)');
            ylabel(obj.axesHandle(2), 'response');
        end

        function handleEpoch(obj, epoch)
            
            % Arrange spikes
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            
            % Count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace); 
            
            % Load frame monitor data
            FMresponse = epoch.getResponse(obj.frameMonitor);
            FMdata = FMresponse.getData();
            times = edu.washington.riekelab.freedland.utils.getFrameTiming(FMdata,0) * 1e-4; % in sec
            
            % Arrange by radial and distance terms
            obj.spins = mod(obj.spinningFrequency .* times .* 360,180); % symmetric, no need to consider separately
            obj.growth = (obj.expandingFrequency .* times) + obj.ovalWidth; % vector of growth
            newEpochResponse = zeros(size(times));
            
            % Downsample to monitor frame rate
            arrangeTimes = [0; times] .* 1e4;
            for a = 1:length(arrangeTimes)-1
                newEpochResponse(a,1) = sum(S.sp >= arrangeTimes(a) & S.sp < arrangeTimes(a+1));
            end
            
            % Now spikes are associated with the spin and growth of the ellipse.
            obj.allAngles = cat(1,obj.allAngles,round(deg2rad(obj.spins),1)); % 5 deg resolution
            obj.allRadii = cat(1,obj.allRadii,round(obj.growth,1)); % 1 px resolution
            obj.allResponses = cat(1,obj.allResponses,newEpochResponse);

            obj.summaryData.allAngles = unique(obj.allAngles);
            obj.summaryData.allRadii = unique(obj.allRadii) .* obj.umPerPix; % convert to um (pix * (um/pix))
            
            % All responses with the same angle
            obj.summaryData.meanAngleResponses = zeros(size(obj.summaryData.allAngles));
            for SpotSizeIndex = 1:length(obj.summaryData.allAngles)
                pullIndices = (obj.summaryData.allAngles(SpotSizeIndex) == obj.allAngles);
                obj.summaryData.meanAngleResponses(SpotSizeIndex) = mean(obj.allResponses(pullIndices));
            end
            obj.lineHandle(1) = polar([obj.summaryData.allAngles; obj.summaryData.allAngles+pi; 2*pi], ...
                [obj.summaryData.meanAngleResponses; obj.summaryData.meanAngleResponses; obj.summaryData.meanAngleResponses(1,1)],'k',...
                'Parent', obj.axesHandle(1));

            % All responses with the same radius
            obj.summaryData.meanRadiiResponses = zeros(size(obj.summaryData.allRadii));
            for SpotSizeIndex = 1:length(obj.summaryData.allRadii)
                pullIndices = (obj.summaryData.allRadii(SpotSizeIndex) == obj.allRadii .* obj.umPerPix); 
                obj.summaryData.meanRadiiResponses(SpotSizeIndex) = mean(obj.allResponses(pullIndices));
            end
            obj.lineHandle(2) = plot(obj.summaryData.allRadii, obj.summaryData.meanRadiiResponses,'Color','black',...
                'Parent', obj.axesHandle(2));
            
        end
        
    end
    
end

