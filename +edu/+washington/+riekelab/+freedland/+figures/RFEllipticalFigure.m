% Adapted from MHT's 'area summation figure'

classdef RFEllipticalFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        preTime
        stimTime
        spotDistances
        resolution
    end
    
    properties (Access = private)
        axesHandle
        location
        distance
        angles
        lineHandle
        fitLineHandle
        allEpochResponses
        baselines
        allAngles
        summaryData
    end
    
    methods
        
        function obj = RFEllipticalFigure(ampDevice, varargin)
            obj.ampDevice = ampDevice;            
            ip = inputParser();
            ip.addParameter('preTime', [], @(x)isvector(x));
            ip.addParameter('stimTime', [], @(x)isvector(x));
            ip.addParameter('spotDistances', 1);
            ip.addParameter('resolution', 90);
            ip.parse(varargin{:});
            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            obj.spotDistances = ip.Results.spotDistances;
            obj.resolution = ip.Results.resolution;
            
            obj.createUi();
        end
        
        function createUi(obj)
            import appbox.*;

            obj.axesHandle = axes( ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle, 'orientation');
            ylabel(obj.axesHandle, 'response');
            obj.setTitle('orientation');
        end
        
        function setTitle(obj, t)
            set(obj.figureHandle, 'Name', t);
            title(obj.axesHandle, t);
        end

        function handleEpoch(obj, epoch)
            % arrange spikes
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            prePts = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            epochResponseTrace = epochResponseTrace(1:prePts+stimPts);
            
            % count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            newEpochResponse = sum(S.sp > prePts); % spike count during stim
            
            a = 0:obj.resolution:360-(1E-9); % angles
            obj.angles = round(deg2rad(a),3);
            obj.location = epoch.parameters('angle'); % individual angle
            obj.distance = epoch.parameters('distance'); % individual distance
            
            obj.allAngles = cat(1,obj.allAngles,round(deg2rad(obj.location),3),round(deg2rad(obj.location)+pi,3)); % add to other size
            obj.allEpochResponses = cat(1,obj.allEpochResponses,newEpochResponse,newEpochResponse);
            
            obj.summaryData.allAngles = unique(obj.allAngles);
            obj.summaryData.meanResponses = zeros(size(obj.summaryData.allAngles));
            for SpotSizeIndex = 1:length(obj.summaryData.allAngles)
                pullIndices = (obj.summaryData.allAngles(SpotSizeIndex) == obj.allAngles);               
                obj.summaryData.meanResponses(SpotSizeIndex) = mean(obj.allEpochResponses(pullIndices));
            end
            
            obj.lineHandle = polar(obj.summaryData.allAngles, obj.summaryData.meanResponses,'--ko',...
                'Parent', obj.axesHandle);
        end
        
    end
    
end

