% Flashes various luminances in center and surround with spatial contrast
% in center
classdef spatialCenterSurroundFlashes < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime     = 250 % ms
        stimTime    = 250 % ms
        tailTime    = 250 % ms
        
        % Center-surround information
        centerRadius            = 100; % radius of receptive-field center, in um
        annulusSize             = 50;  % null region between center and surround, in um
        
        % Luminance information
        centerIntensity         = [-20:10:40];       % luminance of uniform center, in percent contrast
        surroundIntensity       = [-90 -50 0 50 90]; % luminance of uniform surround, in percent contrast
        backgroundIntensity     = 0.168;             % light intensity between trials (0-1)
        
        % Spatial contrast
        centerSlices            = 4; % number of regions to divide center. set to 1 to ignore.
        centerSpatialContrast   = [0 10 20 30 40]; % percent contrast of spatial features in center (+/- % of background)
        
        % Additional settings
        randomizeOrder  = true
        onlineAnalysis  = 'none'
        numberOfAverages = uint16(3) % number of repeats
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        imageTracker
        infoTracker
        counter
        order
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end
        
        function prepareRun(obj)
            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);

            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Generate all spatial masks
            masks = generateMasks(obj);
            maskTracker = zeros(size(masks{1,1}));
            for a = 1:size(masks,1)
                maskTracker = maskTracker + masks{a,1};
            end
            
            % Convert percent contrast to luminance values
            centerLuminances = obj.backgroundIntensity .* (obj.centerIntensity/100 + 1);
            spatialContrasts = obj.backgroundIntensity .* (obj.centerSpatialContrast/100);
            surroundLuminances = obj.backgroundIntensity .* (obj.surroundIntensity/100 + 1);
            
            % Generate images
            obj.imageTracker = cell(length(centerLuminances).*length(spatialContrasts).*...
                length(surroundLuminances),1);
            obj.infoTracker = cell(length(centerLuminances).*length(spatialContrasts).*...
                length(surroundLuminances),3);
            countr = 1;
            for a = 1:length(centerLuminances)
                for b = 1:length(spatialContrasts)
                    for c = 1:length(surroundLuminances)
                        
                        % For tracking parameters
                        obj.infoTracker{countr,1} = obj.centerIntensity(a);
                        obj.infoTracker{countr,2} = obj.centerSpatialContrast(b);
                        obj.infoTracker{countr,3} = obj.surroundIntensity(c);

                        % Generate specific image
                        tmpImage = zeros(size(maskTracker));
                        for d = 1:obj.centerSlices % Center slices
                            if mod(d,2) == 1
                                tmpImage = tmpImage + masks{d,1} .* (centerLuminances(a) + spatialContrasts(b)); % + contrast
                            else
                                tmpImage = tmpImage + masks{d,1} .* (centerLuminances(a) - spatialContrasts(b)); % - contrast
                            end
                        end
                        tmpImage = tmpImage + masks{end,1} .* surroundLuminances(c); % Set surround
                        tmpImage(~maskTracker) = obj.backgroundIntensity;            % Set rest to background intensity
                        
                        obj.imageTracker{countr,1} = uint8(tmpImage .* 255);
                        countr = countr + 1;
                    end
                end
            end
            
            % Setup display
            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(size(obj.imageDatabase,4));
            else
                obj.order = 1:size(obj.imageDatabase,4);
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            % Save parameters as metadata
            epoch.addParameter('centerIntensity',obj.infoTracker{obj.order(obj.counter+1),1});
            epoch.addParameter('centerSpatialContrast',obj.infoTracker{obj.order(obj.counter+1),2});
            epoch.addParameter('surroundIntensity',obj.infoTracker{obj.order(obj.counter+1),3});
        end

        function p = createPresentation(obj)

            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();     
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % Rotate image
            specificImage = obj.imageTracker{obj.order(obj.counter+1),1};
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(specificImage);
            scene.size = [canvasSize(1),canvasSize(2)];
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            % Only allow image to be visible during specific time
            p.addStimulus(scene);
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);

            obj.counter = mod(obj.counter + 1,length(obj.order));
        end    
        
        function masks = generateMasks(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            centerRadiusPx = obj.rig.getDevice('Stage').um2pix(obj.centerRadius); % Convert from um to pixels
            annulusRadiusPx = obj.rig.getDevice('Stage').um2pix(obj.annulusSize); % Convert from um to pixels
            masks = cell(obj.centerSlices+1,1);
            
            % Define space as polar coordinates (r = radial, th = theta)
            [xx,yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            r = sqrt((xx - canvasSize(1)/2).^2 + (yy - canvasSize(2)/2).^2); 
            th = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2));
            th = abs(th-pi/2);              

            % Adjust theta space for strange monitors
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);
            
            % Generate masks
            theta = 0:360/obj.centerSlices:360;
            for a = 1:length(theta)-1
                masks{a,1} = r < (centerRadiusPx-annulusRadiusPx/2) .* (th >= theta(a) & th < theta(a+1));
            end
            masks{end,1} = r > (centerRadiusPx+annulusRadiusPx/2) & r < max(canvasSize)/2;
        end
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages * length(obj.centerSequence));
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages * length(obj.centerSequence));
        end
        
    end
    
end