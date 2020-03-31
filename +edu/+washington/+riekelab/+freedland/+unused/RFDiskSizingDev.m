% Uses dampening to find good disk sizes for retinal projections and
% metamers.
classdef RFDiskSizingDev < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime = 250 % ms
        stimTime = 250 % ms
        tailTime = 250 % ms
        backgroundIntensity = 0.168; % background light intensity
        absoluteMaximum = 0.6;
        
        % RF field information
        rfSigmaCenter = 70; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 170; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Flashes
        cellClass = 'ON'
        contrastSensitivity = 0.2; % contrast for each step
        steps = 3; % number of excitatory steps to make.
        diskWidths = 10:10:150;
        
        % Options
        randomizeOrder = false
        onlineAnalysis = 'none'
        numberOfAverages = uint16(2) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON', 'OFF'})  
        lightSteps
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end
        
        function prepareRun(obj)
            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);

            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.RFDiskSizingFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Pull variables
            obj.micronsPerPixel = obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.monitorSize = obj.rig.getDevice('Stage').getCanvasSize(); % Calculate screen size
            obj.monitorSize = fliplr(obj.monitorSize); % Adjust to [height, width]
            obj.videoSize = floor(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.monitorSize,obj.micronsPerPixel,'PIX2VH')) + 1;
            obj.radii = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.diskWidths,obj.micronsPerPixel,'UM2VH'));

            % Calculate raw intensities
            monitorSteps = (1:obj.steps) .* obj.contrastSensitivity;
            if strcmp(obj.cellClass,'ON')
                monitorSteps = (monitorSteps + 1) * obj.backgroundIntensity;
            elseif strcmp(obj.cellClass,'OFF')
                monitorSteps = (1 - monitorSteps) * obj.backgroundIntensity;
            end

            % Oscillate between lowest and highest intensities
            totalSteps = (1:obj.steps);
            totalStepsR = totalSteps(totalSteps > 1 & totalSteps < obj.steps);
            totalSteps = repmat([totalSteps fliplr(totalStepsR)],1,length(obj.diskWidths));

            % Output light steps as linear equivalent region
            obj.lightSteps = monitorSteps(totalSteps(1:length(obj.diskWidths)));

            if sum(obj.lightSteps < 0) > 1
                error('negative light level encountered. please reduce contrast or number of light steps.')
            end

            % Calculate filter
            rfFilter = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.calculateFilter(obj);

            % Calculate radial coordinates
            [xx,yy] = meshgrid(1:obj.videoSize(2),1:obj.videoSize(1));
            r = sqrt((xx - obj.videoSize(2)/2).^2 + (yy - obj.videoSize(1)/2).^2); 
            th = atan((xx - obj.videoSize(2)/2) ./ (yy - obj.videoSize(1)/2));
            th = abs(th-pi/2);

            % Adjust theta space for strange monitors
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);

            % Build radial filter for image generations
            cuts = 30;
            angles = 0:180/cuts:360;

            radialFilt = zeros(obj.videoSize);
            for a = 1:2:length(angles)-1
                radialFilt(th > angles(a) & th < angles(a+1)) = 3;
            end
            radialFilt(radialFilt == 0) = 2;

            % Identify raw means for image generation
            generatedImg = zeros(obj.videoSize);
            generatedVid = zeros(obj.videoSize(1),obj.videoSize(2),1,length(obj.radii));
            generatedVid_Disk = zeros(obj.videoSize(1),obj.videoSize(2),1,length(obj.radii));
            regionTracker = zeros(obj.videoSize);
            norm_LightStep = zeros(length(obj.radii),1);
            actual_LightStep = zeros(length(obj.radii),1);
            for a = 1:length(obj.radii)

                % Radial region we show to the cell
                region = r < obj.radii(a);

                % Convolve region with RF
                regionSize = sum(region(:));
                T = rfFilter .* region;
                R = sum(T(:)) / regionSize;

                % Region to match
                norm_LightStep(a) = obj.lightSteps(a) .* R;

                % Isolate our new region (we can only change new pixels within this region)
                newRegion = region - regionTracker;
                actualRegion = (newRegion .* radialFilt);

                % We want the entire region to have the correct intensity
                newPixels_1 = sum(actualRegion(:) == 3) / sum(region(:)); % percentage of brighter pixels we can change
                newPixels_2 = sum(actualRegion(:) == 2) / sum(region(:)); % percentage of darker pixels we can change
                oldPixels = 1 - (newPixels_1+newPixels_2); % percentage of other pixels we cannot change

                if a > 1
                    oldAverage = norm_LightStep(a-1);
                else
                    oldAverage = obj.backgroundIntensity;
                end

                % We have two weights of pixels we can apply: bright and dark. We
                % deweight from any pixels we cannot change
                deweightedOldPixels = norm_LightStep(a) - (oldAverage * oldPixels);

                % Actual equation: newPixels_1 * X + newPixels_2 * Y = deweightedOldPixels
                %                  Solve in terms of X (brighter pixels, arbitrarily chosen)
                y = 0:0.01:obj.absoluteMaximum;
                x = (deweightedOldPixels/newPixels_1 - newPixels_2/newPixels_1 * y);

                % Find maximally different region
                x1 = find(x > 0 & x < obj.absoluteMaximum);

                if ~isempty(x1)
                    [~,x1] = max(abs(x(x1)-y(x1)));

                    actualRegion(actualRegion == 3) = x(x1);
                    actualRegion(actualRegion == 2) = y(x1);

                else % No solution exists, need more pixels

                    if sum(x < 0) == length(x)
                        actualRegion(actualRegion == 3) = 0;
                        actualRegion(actualRegion == 2) = 0;
                    end

                    if sum(x > 0) == length(x)
                        actualRegion(actualRegion == 3) = obj.absoluteMaximum;
                        actualRegion(actualRegion == 2) = obj.absoluteMaximum;
                    end
                end

                generatedImg = generatedImg + actualRegion;

                % Check entire region's average
                A = generatedImg .* region;
                actual_LightStep(a) = (sum(A(:)) ./ sum(region(:))) ./ R; % Undo normalization (from perspective of cell)

                regionTracker = regionTracker + newRegion;
                generatedVid(:,:,1,a) = A;
                generatedVid_Disk(:,:,1,a) = region .* (sum(A(:)) ./ sum(region(:)) ./ R);
            end
        end
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages*obj.trials);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages*obj.trials);
        end
        
    end
    
end