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
        valuedImage
        location
        distance
        subpl
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
        end

        function handleEpoch(obj, epoch)
            
            % arrange spikes
            response = epoch.getResponse(obj.ampDevice);
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            prePts = sampleRate*obj.preTime/1000;
            stimPts = sampleRate*obj.stimTime/1000;
            epochResponseTrace = epochResponseTrace(1:prePts+stimPts);
            
            %count spikes
            S = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(epochResponseTrace);
            newEpochResponse = sum(S.sp > prePts); %spike count during stim
            
            if isempty(obj.distance) % first run thru
                obj.distance = 0;
                obj.valuedImage = zeros(200,200); % img displayed during analysis
                obj.subpl = 0;
            end
            
            if obj.distance ~= epoch.parameters('distance')
                obj.subpl = obj.subpl + 1; % new subplot
                obj.valuedImage = zeros(200,200);
            end
            
            a = 0:obj.resolution:180;
            a(a==180) = []; % vector of angles
            
            distSteps = size(obj.spotDistances,2); % number of distances
            obj.location = epoch.parameters('angle'); % individual angle
            obj.distance = epoch.parameters('distance'); % individual distance
            
            if obj.subpl > distSteps
                obj.subpl = 1;
            end
            
            % arrange into cartesian coordinates
            offset = 85;
            xx = round(offset*cos(deg2rad(-obj.location)) + (offset+1));
            yy = round(offset*sin(deg2rad(-obj.location)) + (offset+1));
            xx2 = round(offset*cos(deg2rad(-obj.location + 180)) + (offset+1));
            yy2 = round(offset*sin(deg2rad(-obj.location + 180)) + (offset+1));
            
            sqSize = 15; % size of square to display
            
            % we will now build our image.
            imageAnalysis = zeros(200,200); % current image
            imageAnalysis(yy:yy+sqSize,xx:xx+sqSize) = newEpochResponse;
            imageAnalysis(yy2:yy2+sqSize,xx2:xx2+sqSize) = newEpochResponse;
            
            % define colors
            cmap = colormap(jet);
            
            % make image
            obj.valuedImage = obj.valuedImage + imageAnalysis;
            
            subplot(ceil((distSteps + 1)/2),ceil((distSteps+1)/2),obj.subpl)
            imshow(uint8(obj.valuedImage))
            title(obj.distance)
            colormap(cmap)
        end
        
    end
    
end

