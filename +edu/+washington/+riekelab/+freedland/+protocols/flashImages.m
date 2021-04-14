% Flash rotated naturalistic images.
% By J. Freedland, 2019.
classdef flashImages < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime     = 250 % in ms
        stimTime    = 250 % in ms
        tailTime    = 250 % in ms
        
        % Region of image to be shown
        maskRadius  = 1000;  % in um. Places mask over surrounding portion of image.
        region      = 'full-field'; % where to display image
        
        % Natural image information
        imageNo     = [5,9,2,18,100,56,77,42,56,42,2,78,71,100,6,...
                        100,14,89,96,7,70,78,81,2,64,5,56,2,81,5,55,96,71];  % natural image number (1 to 101)
        frame       = [190,426,957,797,593,550,207,274,943,376,124,298,440,779,603,...
                        357,286,222,171,232,518,259,678,556,437,251,943,621,440,532,629,208,356]; % frame # according to DOVES database. (1 to ~1000).
        observerNo  = 1; % observer number (1 to 19).
        backgroundIntensity = 0.168; % common luminance to hold images at.
        
        % Reduce image appropriately
        includeReducedImages    = true;  % randomize each rotation
        rfSigmaCenter           = 50;    % receptive-field center (only needed if reducing images)
        rfSigmaSurround         = 160;  % receptive-field surround (only needed if reducing images)

        % Additional parameters
        randomize   = true;  % randomize order of flashed images
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(3) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        regionType = symphonyui.core.PropertyType('char', 'row', {'center-only', 'surround-only', 'full-field'}) 
        imageDatabase
        tracker
        counter
        order
        imageTracker
        frameTracker
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
            obj.tracker = repmat({'image'},size(imageFrame,4),1);
            obj.imageTracker = obj.imageNo;
            obj.frameTracker = obj.frame;

            if obj.includeReducedImages == true
                % Load generic settings
                settings = edu.washington.riekelab.freedland.videoGeneration.demoUtils.loadSettings(...
                    obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                    obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'),...
                    obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'),...
                    obj.rfSigmaCenter,obj.rfSigmaSurround);
                
                % Specific settings for reducing image
                settings.diskRegions = settings.diskRadii([1 3 5]); % Place one disk between regions [1] and [3]
                settings.meanDisks   = [1 2];                    % Make disk a linear-equivalent projection
                settings.slices      = 8;
                settings.sliceDisks  = [1 2];                    % Apply slicing to our single disk
                settings.backgroundIntensity = obj.backgroundIntensity*255;
                
                RFFilter = edu.washington.riekelab.freedland.videoGeneration.rfUtils.calculateFilter(settings);
                convolvedImage = imageFrame .* repmat(RFFilter,1,1,1,size(imageFrame,4));
                projection = edu.washington.riekelab.freedland.videoGeneration.utils.linearEquivalency(settings, convolvedImage, RFFilter, imageFrame);
                obj.imageDatabase = cat(4,obj.imageDatabase,uint8(projection));
                obj.tracker = cat(1,obj.tracker,repmat({'reduced'},size(projection,4),1));
                obj.imageTracker = cat(2,obj.imageNo,obj.imageNo);
                obj.frameTracker = cat(2,obj.frame,obj.frame);
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
            epoch.addParameter('imageDisplayed',obj.imageTracker(obj.order(obj.counter+1)));
            epoch.addParameter('frameDisplayed',obj.frameTracker(obj.order(obj.counter+1)));
            epoch.addParameter('flashType',obj.tracker{obj.order(obj.counter+1),1});
            
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
        
        function [imageFrame,backgroundIntensities] = findImage(obj)
            % Relevant frame size
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize() / 2; % radius  
            pixelRange = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(canvasSize,...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'pix2arcmin'));

            imageFrame = zeros(pixelRange(2)*2+1,pixelRange(1)*2+1,1,length(obj.imageNo));
            backgroundIntensities = zeros(length(obj.imageNo),1);
            for a = 1:length(obj.imageNo)
                [path, image] = edu.washington.riekelab.freedland.videoGeneration.utils.pathDOVES(obj.imageNo(a), 1);
                image = image./max(image(:));
                backgroundIntensities(a) = mean(image(:));
                image = image.*255;      

                % Pull image
                imageFrame(:,:,1,a) = image(path.y(obj.frame(a))-pixelRange(2):path.y(obj.frame(a))+pixelRange(2),...
                    path.x(obj.frame(a))-pixelRange(1):path.x(obj.frame(a))+pixelRange(1));
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