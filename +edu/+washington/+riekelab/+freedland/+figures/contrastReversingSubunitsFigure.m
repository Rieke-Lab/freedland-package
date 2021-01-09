% Adapted from MHT's 'area summation figure'

classdef contrastReversingSubunitsFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        ampDevice
        preTime
        stimTime
        temporalFrequency
        monitorSampleRate
        centerRadius
        subunitRadius
    end
    
    properties (Access = private)
        axesHandle
        storeTracker
        storeF1
        storeF2
        storeRatio
        dataTracker
    end
    
    methods
        
        function obj = contrastReversingSubunitsFigure(ampDevice, varargin)
            obj.ampDevice = ampDevice;            
            ip = inputParser();
            ip.addParameter('preTime', [], @(x)isvector(x));
            ip.addParameter('stimTime', [], @(x)isvector(x));
            ip.addParameter('temporalFrequency', [], @(x)isvector(x));
            ip.addParameter('monitorSampleRate', [], @(x)isvector(x));
            ip.addParameter('centerRadius', [], @(x)isvector(x));
            ip.addParameter('subunitRadius', [], @(x)isvector(x));
            ip.parse(varargin{:});
            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            obj.temporalFrequency = ip.Results.temporalFrequency;
            obj.monitorSampleRate = ip.Results.monitorSampleRate;
            obj.centerRadius = ip.Results.centerRadius;
            obj.subunitRadius = ip.Results.subunitRadius;
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
            epochResponseTrace = zeros(1,stimPts);
            epochResponseTrace(S.sp) = 1;

            % Convert to F1/F2 frequency (courtesy of M. Turner)
            L = length(epochResponseTrace); %length of signal, datapoints
            X = abs(fft(epochResponseTrace));
            X = X(1:L/2);
            f = obj.monitorSampleRate*(0:L/2-1)/L; %freq - hz
            [~, F1ind] = min(abs(f-obj.temporalFrequency)); %find index of F1 and F2 frequencies
            [~, F2ind] = min(abs(f-2*obj.temporalFrequency));
            F1 = 2*X(F1ind); %pA^2/Hz for current rec, (spikes/sec)^2/Hz for spike rate
            F2 = 2*X(F2ind); %double b/c of symmetry about zero
            
            % Build mask
            currentSpotSize = epoch.parameters('subunitCoordinates_cartesianMicrons');
            xCoord = obj.centerRadius + currentSpotSize(1);
            yCoord = obj.centerRadius + currentSpotSize(2);
            [xx,yy] = meshgrid(1:ceil(obj.centerRadius)*2,1:ceil(obj.centerRadius)*2);
            r_subunit = sqrt((xx - xCoord).^2 + (yy - yCoord).^2) <= obj.subunitRadius;
            ticklabels = (0:20:size(r_subunit,1)) - obj.centerRadius;
            ticks = (0:20:size(r_subunit,1)*2);
            
            % Save to variable for identification
            obj.dataTracker = cat(1,obj.dataTracker,[epoch.parameters('subunitCoordinates_polarMicrons'), F1, F2]);

            % Assign values 
            if isempty(obj.storeTracker)
                obj.storeTracker    = zeros(size(r_subunit));
                obj.storeF1         = zeros(size(r_subunit));
                obj.storeF2         = zeros(size(r_subunit));
                obj.storeRatio      = zeros(size(r_subunit));
            end
            obj.storeTracker    = obj.storeTracker + r_subunit;
            obj.storeF1         = obj.storeF1 + (F1 .* r_subunit);
            obj.storeF2         = obj.storeF2 + (F2 .* r_subunit);
            obj.storeRatio      = obj.storeRatio + (F2/F1) .* r_subunit;
            
            % Plot
            figure(20)
            subplot(1,3,1)
            imagesc(obj.storeF1 ./ obj.storeTracker);
            title('F1')
            xlabel('x-coordinate (um)')
            ylabel('y-coordinate (um)')
            set(gca, 'XTick', ticks, 'XTickLabel', ticklabels) % in um
            set(gca, 'YTick', ticks, 'YTickLabel', ticklabels)

            subplot(1,3,2)
            imagesc(obj.storeF2 ./ obj.storeTracker);
            title('F2')
            set(gca, 'XTick', ticks, 'XTickLabel', ticklabels)
            set(gca, 'YTick', ticks, 'YTickLabel', ticklabels)
            
            subplot(1,3,3)
            imagesc(obj.storeRatio);
            title('F2/F1')
            set(gca, 'XTick', ticks, 'XTickLabel', ticklabels)
            set(gca, 'YTick', ticks, 'YTickLabel', ticklabels)
        end
    end
end

