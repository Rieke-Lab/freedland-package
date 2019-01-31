classdef amplifyFixations < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        imageNo = 1; % natural image number
        observerNo = 1; % observer number
        
        offsetWidth = 0 % in microns
        offsetHeight = 0 % in microns
        
        amplification = [0 1 10]; % sequence of amplifications
        mirroring = true; % mirror image at edges
        
        onlineAnalysis = 'none'
        numberOfAverages = uint16(2) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        backgroundIntensity
        baseMovement
        baseMovementSequence
        fixMovement
        fixMovementSequence
        pictureInformation
        imageMatrix
        xTraj
        yTraj
        timeTraj
        currentStimSet
        randomIndex
        selectionIndex
        xTrajAmp
        yTrajAmp
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)
                   
            % For Symphony
            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Identify image
            imageIdentifier = [5 8 13 17 23 27 31 39 56 64 100];
            imageVal = imageIdentifier(obj.imageNo);
            
            % Grab image information
            [~, ~, ~, obj.pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal, obj.observerNo,...
                 'amplification', 0, 'mirroring', obj.mirroring);
            
            % Pre-allocate memory
            obj.baseMovementSequence.x = zeros(size(obj.amplification,2),size(obj.pictureInformation.saccadeTracking,2));
            obj.baseMovementSequence.y = zeros(size(obj.amplification,2),size(obj.pictureInformation.saccadeTracking,2));
            obj.fixMovementSequence.x = zeros(size(obj.amplification,2),size(obj.pictureInformation.saccadeTracking,2));
            obj.fixMovementSequence.y = zeros(size(obj.amplification,2),size(obj.pictureInformation.saccadeTracking,2));
            
            % Convert um to pixels
            obj.offsetHeight = obj.rig.getDevice('Stage').um2pix(obj.offsetHeight);
            obj.offsetWidth = obj.rig.getDevice('Stage').um2pix(obj.offsetWidth);
            
            % Make paths for every amplification.
            for a = 1:size(obj.amplification,2)
                [~, obj.baseMovement, obj.fixMovement, ~] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal, obj.observerNo,...
                    'amplification', obj.amplification(1,a),'offSetHeight', obj.offsetHeight,'offSetWidth',...
                    obj.offsetWidth, 'mirroring', obj.mirroring);
                obj.baseMovementSequence.x(a,:) = obj.baseMovement.x;
                obj.baseMovementSequence.y(a,:) = obj.baseMovement.y;
                obj.fixMovementSequence.x(a,:) = obj.fixMovement.x;
                obj.fixMovementSequence.y(a,:) = obj.fixMovement.y;
            end
                
            % Adjust image
            img = obj.pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img = img.*255;
            obj.imageMatrix = uint8(img);
            
            % Recombine trajectories
            obj.xTraj = obj.baseMovementSequence.x + obj.fixMovementSequence.x;
            obj.yTraj = obj.baseMovementSequence.y + obj.fixMovementSequence.y;
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % DOVES resolution
           
            % Adjust axes and units for monitor
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            obj.xTraj = obj.xTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj = obj.yTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            
            % To avoid adaptation, we randomize the order
            obj.selectionIndex = 1;
            obj.randomIndex = randperm(size(obj.amplification,2));
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('individualAmplification', obj.amplification(obj.randomIndex(obj.selectionIndex)));
            epoch.addParameter('currentStimSet', obj.currentStimSet);
        end
        
        function p = createPresentation(obj)
            
            % Prep stage for presentation
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(obj.imageMatrix);
            scene.size = [size(obj.imageMatrix,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.imageMatrix,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            
            % Select the random index
            obj.xTrajAmp = obj.xTraj(obj.randomIndex(obj.selectionIndex),:);
            obj.yTrajAmp = obj.yTraj(obj.randomIndex(obj.selectionIndex),:);
            
            % Apply eye trajectories to move image around
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3, p0));
            
            function p = getScenePosition(obj, time, p0)
                if time < 0
                    p = p0;
                elseif time > obj.timeTraj(end) %out of eye trajectory, hang on last frame
                    p(1) = p0(1) + obj.xTrajAmp(end);
                    p(2) = p0(2) + obj.yTrajAmp(end);
                else % within eye trajectory and stim time
                    dx = interp1(obj.timeTraj,obj.xTrajAmp,time);
                    dy = interp1(obj.timeTraj,obj.yTrajAmp,time);
                    p(1) = p0(1) + dx;
                    p(2) = p0(2) + dy;
                end
            end            
            
            p.addStimulus(scene);
            p.addController(scenePosition);

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            if obj.selectionIndex < size(obj.amplification,2)
                obj.selectionIndex = obj.selectionIndex + 1; % choose next index.
            elseif obj.selectionIndex == size(obj.amplification,2) && obj.numEpochsPrepared < obj.numberOfAverages * size(obj.amplification,2)
                obj.selectionIndex = 1;
                obj.randomIndex = randperm(size(obj.amplification,2));
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * size(obj.amplification,2);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * size(obj.amplification,2);
        end

    end
    
end