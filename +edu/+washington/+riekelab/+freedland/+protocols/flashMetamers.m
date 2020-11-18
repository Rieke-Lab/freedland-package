% Flash rotated naturalistic images.
% By J. Freedland, 2019.
classdef flashMetamers < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Cell information
        rfSigmaCenter = 30; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % image info
        slices              = 8; % number of slices
        contrast            = 0.5; % 0 - 1
        randomize           = true;
        
        backgroundIntensity = 0.168; % 0 - 1 

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(3) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        backgroundIntensity
        imageDatabase
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
            obj.showFigure('edu.washington.riekelab.freedland.figures.receptiveFieldFitFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','experimentID');
            
            % Load settings
            retinalMetamers = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.demoUtils.loadSettings(70,170);
            retinalMetamers.slices              = obj.slices;
            retinalMetamers.sliceDisks          = 1;
            retinalMetamers.metamerDisks        = 1;
            retinalMetamers.imageNo             = 0;
            retinalMetamers.diskRegions = obj.diskRadii([1 3]); % Place one disk between regions [1] and [3]
            
            % Generate intensity values
            A = nchoosek(1:obj.slices,round(obj.slices/2)); % 
            diskContrasts = ones(8,size(A,1)) * -obj.contrast;
            for a = 1:size(A,1)
                diskContrasts(A(a,:),a) = obj.contrast;
            end
            stimulus.values = obj.backgroundIntensity*255 + (diskContrasts .* (obj.backgroundIntensity*255)); % ON cell

            %%% Make metamers
            % Pull cell-specific receptive field
            [RFFilter,obj.rfSizing] = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.rfUtils.calculateFilter(obj);

            % Compare stimulus against database
            disp('Pulling metamer library...')
            databaseTraj = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.metamerUtils.pullLibrary(obj);
            weightedDatabaseTraj = databaseTraj .* RFFilter; % Convolve with stimulus

            % Find linear-equivalent regions for cell-specific RF
            disp('Calculating low-dimensional projection...')
            [~,databaseValues,stimulus.masks] =  edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.linearEquivalency(obj, weightedDatabaseTraj, RFFilter);

            % Find best match to designed stimulus
            output =  edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.metamerUtils.findReplacements(obj,stimulus,databaseValues,databaseTraj);
            obj.imageDatabase = squeeze(output.metamer);
            %%%
            
            % Setup display
            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(1:size(obj.imageDatabase,3));
            else
                obj.order = 1:size(obj.imageDatabase,3);
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('experimentID',obj.order(obj.counter+1));
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();     
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % Rotate image
            specificImage = obj.imageDatabase(:,:,obj.order(obj.counter+1));
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(specificImage);
            sz = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(...
                size(specificImage,1),obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'VH2PIX');
            scene.size = [sz sz];
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
            
            % Add mask
            aperatureDiameter = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(...
                obj.maskRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'UM2PIX') .* 2);
            
            if (aperatureDiameter > 0) %% Create aperture
                aperture = stage.builtin.stimuli.Rectangle();
                aperture.position = canvasSize/2;
                aperture.color = obj.backgroundIntensity;
                aperture.size = [max(canvasSize) max(canvasSize)];
                mask = stage.core.Mask.createCircularAperture(aperatureDiameter/max(canvasSize), 1024); %circular aperture
                aperture.setMask(mask);
                p.addStimulus(aperture); %add aperture
            end
            
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.order);
        end
    end
end