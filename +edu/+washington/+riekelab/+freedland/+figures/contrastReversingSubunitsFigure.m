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
        summaryData
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
            
            obj.createUi();
        end
        
        function createUi(obj)
            import appbox.*;

            iconDir = [fileparts(fileparts(mfilename('fullpath'))), '\+utils\+icons\'];
            toolbar = findall(obj.figureHandle, 'Type', 'uitoolbar');
            calculateWeights = uipushtool( ...
                'Parent', toolbar, ...
                'TooltipString', 'Export subunits', ...
                'Separator', 'on', ...
                'ClickedCallback', @obj.onSelectedExportSubunits);
            setIconImage(calculateWeights, [iconDir, 'DoG.png']);
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
            currentSpotSize = epoch.parameters('pixelCoordinates_polarPx');
            xCoord = obj.centerRadius + currentSpotSize(1) .* cos(deg2rad(currentSpotSize(2)));
            yCoord = obj.centerRadius - currentSpotSize(1) .* sin(deg2rad(currentSpotSize(2)));
            [xx,yy] = meshgrid(1:obj.centerRadius*2,1:obj.centerRadius*2);
            r_subunit = sqrt((xx - xCoord).^2 + (yy - yCoord).^2) <= obj.subunitRadius;
            
            % Save to variable for identification
            obj.dataTracker = cat(1,obj.dataTracker,[epoch.parameters('pixelCoordinates_polarMicrons'), F1, F2]);
                
            % First run
            if ~isfield(obj.summaryData,'F1')
                obj.summaryData.F1 = zeros(obj.centerRadius*2,obj.centerRadius*2);
                obj.summaryData.F2 = zeros(obj.centerRadius*2,obj.centerRadius*2);
                obj.summaryData.ratio = zeros(obj.centerRadius*2,obj.centerRadius*2);
                obj.summaryData.tracker = zeros(obj.centerRadius*2,obj.centerRadius*2);
            end
            
            % Assign values
            obj.summaryData.tracker = obj.summaryData.tracker + r_subunit;
            obj.summaryData.F1 = obj.summaryData.F1 + (F1 .* r_subunit);
            obj.summaryData.F2 = obj.summaryData.F2 + (F2 .* r_subunit);
            obj.summaryData.ratio = obj.summaryData.ratio + (F2/F1) .* r_subunit;

            figure(1)
            subplot(1,3,1)
            imagesc(obj.summaryData.F1 ./ obj.summaryData.tracker);
            title('F1')

            subplot(1,3,2)
            imagesc(obj.summaryData.F2 ./ obj.summaryData.tracker);
            title('F2')
            
            subplot(1,3,3)
            imagesc(obj.summaryData.ratio);
            title('F2/F1')
        end
        
        function onSelectedExportSubunits(obj, ~, ~)
            % Sort by F2 frequency
            [~,i] = sort(obj.dataTracker(:,4),'descend');
            
            % Export as .txt.
            dlmwrite(strcat('Documents/subunitLocation_um.txt'),obj.dataTracker(i,:))
        end
    end
end

