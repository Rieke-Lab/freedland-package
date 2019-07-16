% Flash natural images based on a series of manually-defined presets.
% By J. Freedland, 2019.
classdef RFFlashImages < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        imageNo = 5; % natural image number (1 to 101), as a vector.
        observerNo = 1; % observer number (1 to 19). Can be a single value or vector.
        frameNumber = 200; % frame number according to DOVES database. Must be a vector.
        presetVal = '1'; % can override image number and frame number using presets (if needed).
        randomlyDisplay = true; % display images in a random order

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(10) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        meanIntegrationType = symphonyui.core.PropertyType('char', 'row', {'uniform','gaussian'})
        overrideCoordinateType = symphonyui.core.PropertyType('char', 'row', {'pixels','RF'})
        replacementImageType = symphonyui.core.PropertyType('char', 'row', {'disk array','null array'})
        presetValType = symphonyui.core.PropertyType('char', 'row', {'none','1','2'})
        backgroundIntensity
        imageDatabase
        counter
        permut
        observerVector
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
            
            if strcmp(obj.presetVal,'1')
                obj.imageNo = [5 2 2 12 6 41 66 97 101 5 39 70 97];
                obj.frameNumber = [200 969 148 664 180 854 961 697 144 600 288 96 9];
            elseif strcmp(obj.presetVal,'2')
                obj.imageNo = [3,71,43,63,87,61,60,19,42,33,78,75,77,46,28,36,39,67,68,82,74,53,34,101,9,40,37,49,99,31,24,65,97,54,96,56,66];
                obj.frameNumber = [191,92,545,229,131,530,168,103,793,89,381,365,907,440,771,377,603,180,851,642,96,607,419,6,993,324,116,60,333,273,986,204,828,84,168,803,518];
            end
            
            % Rearrange observers into a vector for ease.
            if length(obj.observerNo) == 1
                obj.observerVector = repelem(obj.observerNo,length(obj.imageNo));
            else
                obj.observerVector = obj.observerNo;
            end
            
            % Catch common errors.
            if length(obj.imageNo) ~= length(obj.frameNumber)
                error('The number of images and frames must be equivalent')
            end
            
            if length(obj.imageNo) ~= length(obj.observerVector)
                error('Incorrect number of observers')
            end
            
            % Use DOVES database to identify images along trajectory
            imageFrames = findImage(obj);
            obj.imageDatabase = uint8(imageFrames);
            
            obj.counter = 0;
            if obj.randomlyDisplay == true
                obj.permut = randperm(length(obj.imageNo));
            else
                obj.permut = 1:length(obj.imageNo);
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            % find appropriate stimulus
            A = obj.permut(obj.counter+1);
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('imageNo', obj.imageNo(A));
            epoch.addParameter('observerNo', obj.observerVector(A));
            epoch.addParameter('frameNo', obj.frameNumber(A));
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();             
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);

            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            A = obj.permut(obj.counter+1);
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(obj.imageDatabase(:,:,A));
            scene.size = [size(obj.imageDatabase,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.imageDatabase,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            p.addStimulus(scene);
            
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            obj.counter = mod(obj.counter + 1,length(obj.imageNo));
            
            if obj.counter == 0 && obj.randomlyDisplay == true
                obj.permut = randperm(length(obj.imageNo)); % randomly rearrange vector at end
            end
        end
        
        % Apply RF Filter over the entire trajectory.
        function imageFrames = findImage(obj)
            
            % Calculate frame size
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            imgSize = ceil(canvasSize / (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')));
            xRange = floor(imgSize(1) / 2);
            yRange = floor(imgSize(2) / 2);
            
            imageFrames = zeros(imgSize(2),imgSize(1),length(obj.imageNo)); % matrix with images
            
            for a = 1:length(obj.imageNo) % Walk along our vectors
                tempImage = obj.imageNo(a);
                tempObserver = obj.observerVector(a);
                [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(tempImage, tempObserver,...
                        'amplification', 1,'mirroring', true); % Pull coordinates of DOVES database
                
                % Scale pixels in image to monitor
                img = pictureInformation.image;
                img = (img./max(max(img)));              
                img2 = img.*255;
                
                if a == 1
                    obj.backgroundIntensity = mean(img(:)); % Set background intensity to reference image
                end
            
                % Produce trajectories
                xTraj = baseMovement.x(obj.frameNumber(a)) + fixMovement.x(obj.frameNumber(a));
                yTraj = baseMovement.y(obj.frameNumber(a)) + fixMovement.y(obj.frameNumber(a));
                
                % While the rig automatically centers the stimulus, our
                % calculation doesn't.
                centering = obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'); % in mu
                centeringPix = centering ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
                centeredXTraj = round(xTraj - centeringPix(1));
                centeredYTraj = round(yTraj + centeringPix(2));
                
                % create image
                imageFrames(:,:,a) = img2(centeredYTraj-yRange:centeredYTraj+yRange-1,...
                    centeredXTraj-xRange:centeredXTraj+xRange-1); 
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.imageNo);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.imageNo);
        end
    end
end