% Flash rotated naturalistic images.
% By J. Freedland, 2019.
classdef flashImages < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Region of image to be shown
        maskRadius = 100;  % in um. Places mask over surrounding portion of image.
        region = 'full-field'; % where to display image
        
        % Natural image
        imageNo = [5,5,12,12,71,71,73,73,79,79,81,81,100,100];              % natural image number (1 to 101)
        frame = [532,471,458,608,900,714,681,946,252,513,678,64,593,779];   % frame # according to DOVES database. (1 to ~1000).
        observerNo = 1;                                                     % observer number (1 to 19).
        randomize = true;  % randomize each rotation

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(3) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        regionType = symphonyui.core.PropertyType('char', 'row', {'center-only', 'surround-only', 'full-field'}) 
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
            
            % Check for errors
            if length(obj.imageNo) ~= length(obj.frame)
                error('The number of images and frames to probe must be equal.')
            end
            
            % Use DOVES database to identify images along trajectory
            imageFrame = findImage(obj);
            obj.imageDatabase = uint8(imageFrame);

            % Setup display
            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(length(obj.imageNo));
            else
                obj.order = 1:length(obj.imageNo);
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('imageDisplayed',obj.imageNo(obj.order(obj.counter+1)));
            epoch.addParameter('frameDisplayed',obj.frame(obj.order(obj.counter+1)));
            
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
            specificImage = obj.imageDatabase(:,:,1,obj.order(obj.counter+1));
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(specificImage);
            sz = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                size(specificImage,1),obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
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
            
            % Create aperature
            aperatureDiameter = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.maskRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix') .* 2);
            aperture = stage.builtin.stimuli.Rectangle();
            aperture.position = canvasSize/2;
            aperture.color = obj.backgroundIntensity;
            aperture.size = [max(canvasSize) max(canvasSize)];
            
            if strcmp(obj.region,'center-only') %% Create aperture
                mask = stage.core.Mask.createCircularAperture(aperatureDiameter/max(canvasSize), 1024);
            elseif strcmp(obj.region,'surround-only')
                mask = stage.core.Mask.createAnnulus(0,aperatureDiameter/max(canvasSize), 1024);
            elseif strcmp(obj.region,'full-field')
                mask = stage.core.Mask.createCircularAperture(1, 1024);
            end
            
            aperture.setMask(mask);
            p.addStimulus(aperture);
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end
        
        function imageFrame = findImage(obj)
            
            % Relevant frame size
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();   
            pixelRange = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(max(canvasSize),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'pix2arcmin'));

            imageFrame = zeros(pixelRange*2+1,pixelRange*2+1,1,length(obj.imageNo));
            for a = 1:length(obj.imageNo)
                [path, image] = edu.washington.riekelab.freedland.videoGeneration.utils.pathDOVES(obj.imageNo(a), 1);
                image = image./max(image(:));
                obj.backgroundIntensity = mean(image(:));
                image = image.*255;      

                % Pull image
                imageFrame(:,:,1,a) = image(path.y(obj.frame(a))-pixelRange:path.y(obj.frame(a))+pixelRange,...
                    path.x(obj.frame(a))-pixelRange:path.x(obj.frame(a))+pixelRange);
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.order);
        end
    end
end