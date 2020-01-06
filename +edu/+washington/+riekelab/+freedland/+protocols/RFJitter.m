% Compares mirrored naturalistic image trajectories with and without
% fixational eye movements.
% By J. Freedland, 2019.
classdef RFJitter < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        imageNo = 5; % natural image number (1 to 101), as a vector.
        observerNo = 1; % observer number (1 to 19). Can be a single value or vector.
        fixationIntensity = [0 1 2];
        randomize = true; % display images in a random order

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
        
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        backgroundIntensity
        sequence
        xTraj
        yTraj
        flippedxTraj
        imageMatrix
        flippedImageMatrix
        selection
        timeTraj
        counter
        doubleCounter
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
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',2);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Calculate basic trajectory
            [~, baseMovement, ~, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.imageNo, obj.observerNo,...
                        'amplification', 1,'mirroring', true);
                    
            % Scale image pixels to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix = uint8(img2);
                    
            obj.xTraj = zeros(length(baseMovement.x),length(obj.fixationIntensity));
            obj.yTraj = zeros(length(baseMovement.x),length(obj.fixationIntensity));
            for a = 1:length(obj.fixationIntensity)
                
                [~, ~, fixMovement, ~] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.imageNo, obj.observerNo,...
                        'amplification', obj.fixationIntensity(a),'mirroring', true);
                    
                obj.xTraj(:,a) = baseMovement.x + fixMovement.x;
                obj.yTraj(:,a) = baseMovement.y + fixMovement.y;
            end

            % We do not need to consider the entire trajectory, however.
            frames = round((obj.stimTime + 50) / 1000 * 200); % max # of frames in DOVES database, with 50ms cushion

            if size(obj.xTraj,1) <= frames
                frames = size(obj.xTraj,1);
            end
            
            obj.xTraj = obj.xTraj(1:frames,:);
            obj.yTraj = obj.yTraj(1:frames,:);
            obj.flippedxTraj = abs(obj.xTraj - size(obj.imageMatrix,2)); % Mirror x coordinates
            obj.timeTraj = (0:(size(obj.xTraj,1)-1)) ./ 200; % DOVES resolution
                
            % Scale image pixels to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix = uint8(img2);
            obj.flippedImageMatrix = uint8(fliplr(img2)); % Flip image for reference
                
            % Invert for monitor
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            obj.flippedxTraj = -(obj.flippedxTraj - size(img,2)/2);
            
            % Convert from DOVES (VH units) to monitor (pixel units).
            obj.xTraj = obj.xTraj .* (3.3 ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            obj.yTraj = obj.yTraj .* (3.3 ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            obj.flippedxTraj = obj.flippedxTraj .* (3.3 ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));

            if obj.randomize == false
                obj.sequence = 1:length(obj.fixationIntensity);
            else
                obj.sequence = randperm(length(obj.fixationIntensity));
            end
            obj.counter = 1;
            obj.doubleCounter = 1;
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime)*2 / 1e3;
            
            % Find appropriate stimulus
            A = obj.sequence(obj.counter);
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('specificFixationIntensity', obj.fixationIntensity(A));
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % in normal pixels            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 2 * 1e-3);
            
            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity)
            
            % Find corresponding experiment
            obj.selection = obj.sequence(obj.counter);
            
            % Base image
            scene = stage.builtin.stimuli.Image(obj.imageMatrix);
            sceneMirrored = stage.builtin.stimuli.Image(obj.flippedImageMatrix);
            
            scene.size = [size(obj.imageMatrix,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.imageMatrix,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            sceneMirrored.size = scene.size;
            
            p0 = canvasSize/2;
            scene.position = p0;
            sceneMirrored.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            sceneMirrored.setMinFunction(GL.LINEAR);
            sceneMirrored.setMagFunction(GL.LINEAR);

            % Program each scene
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3, p0, obj.xTraj(:,obj.selection)));
            mirroredPosition = stage.builtin.controllers.PropertyController(sceneMirrored,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3 - cycleTime, p0, obj.flippedxTraj(:,obj.selection)));

            function p = getScenePosition(obj, time, p0, xTraj)
                
                yTraject = obj.yTraj(:,obj.selection);
                
                if time <= 0
                    p = p0;
                elseif time > obj.timeTraj(end) % Beyond eye trajectory, hang on last frame
                    p(1) = p0(1) + xTraj(end);
                    p(2) = p0(2) + yTraject(end);
                else % Within eye trajectory and stim time
                    dx = interp1(obj.timeTraj,xTraj,time);
                    dy = interp1(obj.timeTraj,yTraject,time);
                    p(1) = p0(1) + dx;
                    p(2) = p0(2) + dy;
                end
            end
            
            p.addStimulus(scene);
            p.addController(scenePosition);
            
            p.addStimulus(sceneMirrored);
            p.addController(mirroredPosition);
            
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            sceneMirroredVisible = stage.builtin.controllers.PropertyController(sceneMirrored, 'visible', ...
                    @(state)state.time >= (obj.preTime) * 1e-3 + cycleTime && state.time < (obj.preTime + obj.stimTime) * 1e-3 + cycleTime);
            
            p.addController(sceneVisible);
            p.addController(sceneMirroredVisible);
            
            if (obj.doubleCounter - 1) == obj.numberOfAverages
                obj.doubleCounter = 1;
                obj.counter = obj.counter + 1;
            else
                obj.doubleCounter = obj.doubleCounter + 1;
            end
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.fixationIntensity);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.fixationIntensity);
        end
    end
end