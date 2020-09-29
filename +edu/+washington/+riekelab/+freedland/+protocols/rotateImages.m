% Flash natural images based on a series of manually-defined presets.
% Added updating background intensity: 08/26/2019
% By J. Freedland, 2019.
classdef rotateImages < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Natural image
        imageNo = 81; % natural image number (1 to 101), as a vector.
        observerNo = 1; % observer number (1 to 19). Can be a single value or vector.
        frame = 200; % frame # according to DOVES database. Typically between 1 and 1000.
        
        % Rotation characteristics
        stepRotation = 30; % degrees to rotate each cycle
        maskRadius = 300; % in um. Places mask over surrounding portion of image.
        randomize = false; % randomize each rotation

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        backgroundIntensity
        imageDatabase
        counter
        order
        rotations
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
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','rotation (deg)');
            
            % Use DOVES database to identify images along trajectory
            [imageFrame, obj.backgroundIntensities] = findImage(obj);
            obj.imageDatabase = uint8(imageFrame);
            
            % Identify rotations
            obj.rotations = unique(mod(0:obj.stepRotation:360,360));
            
            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(length(obj.rotations));
            else
                obj.order = 1:length(obj.rotations);
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('rotation',obj.rotations(obj.order(obj.counter+1)));
            
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
            specificImage = imrotate(obj.imageDatabase,obj.rotations(obj.order(obj.counter+1)));
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
        
        function imageFrame = findImage(obj)
            
            % Relevant frame size
            pixelRange = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(obj.maskRadius,...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'UM2VH'));

            [empiricalPath, ~, ~, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(...
                obj.imageNo,obj.observerNo);
                
            % Scale pixels in image to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));  
            obj.backgroundIntensity = mean(img(:));
            img = img.*255;      
                
            % Pull image
            imageFrame = img(empiricalPath.y(obj.frame)-pixelRange:empiricalPath.y(obj.frame)+pixelRange,...
                empiricalPath.x(obj.frame)-pixelRange:empiricalPath.x(obj.frame)+pixelRange);
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.imageNo);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.imageNo);
        end
    end
end