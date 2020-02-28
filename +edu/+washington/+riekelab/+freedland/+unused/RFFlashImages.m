% Flash natural images based on a series of manually-defined presets.
% Added updating background intensity: 08/26/2019
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
        presetVal = '2'; % can override image number and frame number using presets (if needed).
        randomlyDisplay = true; % display images in a random order

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        presetValType = symphonyui.core.PropertyType('char', 'row', {'none','1','2'})
        backgroundIntensity
        imageDatabase
        counter
        permut
        observerVector
        backgroundIntensities
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
                obj.imageNo = [56,61,5,19,58,30,89,28,44,2,85,60,89,46,95,75,18,35,46,83,79,72,12,31,51,81,86,100,80,22,21,38,85,8,91,32,89,16,20,55,83,42,17,11,93,23,79,5];
                obj.frameNumber = [859,507,331,272,838,506,366,593,1026,437,862,634,583,116,468,90,112,143,831,217,636,311,1040,769,938,594,1012,1013,823,646,187,889,766,894,602,682,594,611,956,932,148,238,880,514,249,91,538,472];
            end
            
            % Rearrange observers into a vector for ease.
            if length(obj.observerNo) == 1
                obj.observerVector = repelem(obj.observerNo,length(obj.imageNo));
            else
                obj.observerVector = obj.observerNo;
            end
            
            % Catch common errors.
            if length(obj.imageNo) ~= length(obj.frameNumber)
                error('Each image must have a frame to show.')
            end
            
            if length(obj.imageNo) ~= length(obj.observerVector)
                error('Incorrect number of observers')
            end
            
            % Use DOVES database to identify images along trajectory
            [imageFrames, obj.backgroundIntensities] = findImage(obj);
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
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % in normal pixels            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);

            A = obj.permut(obj.counter+1);
            
            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensities(A,1))
            
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
        
        function [imageFrames, backgroundIntensities] = findImage(obj)
            
            % Calculate frame size
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            imgSize = ceil(canvasSize / (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'))); % Convert to DOVES VH units
            xRange = floor(imgSize(1) / 2);
            yRange = floor(imgSize(2) / 2);
            
            imageFrames = zeros(yRange.*2+1,xRange.*2+1,length(obj.imageNo)); % Collection of images
            backgroundIntensities = zeros(length(obj.imageNo),1);
            
            for a = 1:length(obj.imageNo)
                tempImage = obj.imageNo(a); % Pull image #
                tempObserver = obj.observerVector(a); % Pull observer #
                [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(tempImage, tempObserver,...
                        'amplification', 1,'mirroring', true); % Pull coordinates from DOVES database
                    
                % Pull frame #
                xTraj = baseMovement.x(obj.frameNumber(a)) + fixMovement.x(obj.frameNumber(a));
                yTraj = baseMovement.y(obj.frameNumber(a)) + fixMovement.y(obj.frameNumber(a));
                
                % Scale pixels in image to monitor
                img = pictureInformation.image;
                img = (img./max(max(img)));              
                img2 = img.*255;      
                
                % Pull image
                imageFrames(:,:,a) = img2(yTraj-yRange:yTraj+yRange,...
                    xTraj-xRange:xTraj+xRange); 
                
                % Identify global average
                R = imageFrames(:,:,a);
                backgroundIntensities(a,1) = mean(R(:)) ./ 255;
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