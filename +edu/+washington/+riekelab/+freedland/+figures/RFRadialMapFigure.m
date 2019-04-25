% Adapted from MHT's 'area summation figure'

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

            obj.axesHandle(1) = subplot(1,3,1, ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle(1), 'orientation');
            ylabel(obj.axesHandle(1), 'response');
            
            obj.axesHandle(2) = subplot(1,3,2, ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle(2), 'distance from center (um)');
            ylabel(obj.axesHandle(2), 'response');
            
            obj.axesHandle(3) = subplot(1,3,3, ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle(3), 'pixel');
            ylabel(obj.axesHandle(3), 'pixel');
        end

        function handleEpoch(obj, epoch)
            
            % arrange spikes
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            
            % count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace); 
            
            %load frame monitor data
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
            
            % Now the times are associated with the spin and growth of the
            % ellipse.
            obj.allAngles = cat(1,obj.allAngles,round(deg2rad(obj.spins),2)); % 5 deg resolution
            obj.allRadii = cat(1,obj.allRadii,round(obj.growth,-1)); % 10 px resolution
            obj.allResponses = cat(1,obj.allResponses,newEpochResponse);

            obj.summaryData.allAngles = unique(obj.allAngles);
            obj.summaryData.allRadii = unique(obj.allRadii);
            
            obj.summaryData.meanAngleResponses = zeros(size(obj.summaryData.allAngles));
            for SpotSizeIndex = 1:length(obj.summaryData.allAngles)
                pullIndices = (obj.summaryData.allAngles(SpotSizeIndex) == obj.allAngles);
                
                % all responses with the same angle
                obj.summaryData.meanAngleResponses(SpotSizeIndex) = mean(obj.allResponses(pullIndices));
            end
            
            obj.summaryData.meanRadiiResponses = zeros(size(obj.summaryData.allRadii));
            for SpotSizeIndex = 1:length(obj.summaryData.allRadii)
                pullIndices = (obj.summaryData.allRadii(SpotSizeIndex) == obj.allRadii); 
                
                % all responses with the same radius
                obj.summaryData.meanRadiiResponses(SpotSizeIndex) = mean(obj.allResponses(pullIndices));
            end
            
            obj.lineHandle(1) = polar([obj.summaryData.allAngles; obj.summaryData.allAngles+pi], ...
                [obj.summaryData.meanAngleResponses; obj.summaryData.meanAngleResponses],'k',...
                'Parent', obj.axesHandle(1));
            
            obj.summaryData.allRadii = obj.summaryData.allRadii .* obj.umPerPix; % convert to um (pix * (um/pix))

            obj.lineHandle(2) = line(obj.summaryData.allRadii, obj.summaryData.meanRadiiResponses,'Color','black',...
                'Parent', obj.axesHandle(2));
            
            % plot RF
            xx = zeros(size(obj.spins,1),1);
            yy = zeros(size(obj.spins,1),1);
            for a = 1:size(obj.spins,1)
                xx(a,1) = round(cos(obj.spins(a,1)) .* obj.growth(a,1));
                yy(a,1) = round(sin(obj.spins(a,1)) .* obj.growth(a,1));
            end

            % move radius = 0 to center of image
            xx = xx + round(obj.canvasSize(2)*obj.umPerPix/2);
            yy = yy + round(obj.canvasSize(1)*obj.umPerPix/2);

            zz = diff([0; newEpochResponse]); % change in spikes gives RF
            zz = zz ./ max(zz(:)); % normalize
            
            minXY = floor(mean(sqrt(diff(sort(xx)).^2 + diff(sort(yy)).^2))); % avg pixel size
            
            meshGrid = zeros(round(obj.canvasSize(2)*obj.umPerPix),round(obj.canvasSize(1)*obj.umPerPix)).*0.1;
            for a = 1:size(xx,1)
                meshGrid(xx(a,1):xx(a,1)+minXY,yy(a,1):yy(a,1)+minXY) = zz(a,1);
            end
            
            obj.totalMeshGrid = cat(3,obj.totalMeshGrid,meshGrid);
            obj.totalMeshGrid = mean(obj.totalMeshGrid,3);
            
            obj.lineHandle(3) = imshow(obj.totalMeshGrid,...
                'Parent', obj.axesHandle(3));
            
        end
        
    end
    
end

