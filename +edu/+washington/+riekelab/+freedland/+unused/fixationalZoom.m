classdef fixationalZoom < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 50 % in ms
        stimTime = 200 % in ms
        tailTime = 50 % in ms
        controlTime = 100 % in ms
        
        imageStyle = 'dotsW'; % image to move
        observerNo = 1; % observer number
      
        amplification = [0 0.5 1 2 4 8]; % sequence of amplifications
        zooming = [0.8 0.9 1 2 4 8]; % factor of zoom. <1 is zooming out, >1 zooms in.
        
        mirroring = true; % mirror image on edges
        randomFixations = false; % change DOVES fixation path between each average (full cycle)
        
        onlineAnalysis = 'none'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        imageStyleType = symphonyui.core.PropertyType('char', 'row', {'dotsW','dotsSW','dotsS','dotsSE',...
            'a_dotsW','a_dotsSW','a_dotsS','a_dotsSE'})
        backgroundIntensity
        pictureInformation
        fixMovement
        tempImage
        imageMatrix
        picture
        xTraj
        yTraj
        timeTraj
        currentStimSet
        zoomAmpPerm
        randomIndex
        selectionIndex
        xTrajAmp
        yTrajAmp
        saccadeIndex
        fixationIndexAdder
        fixationIndex
        randomFixationIndex
        imageMatrixFull
        tempImage2
        mirroredImg
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
                        
            % Choose a random DOVES trajectory
            [~, ~, obj.fixMovement, obj.pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(5, obj.observerNo,...
                 'amplification', 1);
             
            % Catalog all fixational eye movements in trajectory
            obj.saccadeIndex = find(obj.pictureInformation.saccadeTracking); % find fixations
            obj.fixationIndexAdder = 1; % for grouping
            
            % Find DOVES fixational eye movements of suitable length.
            for a = 1:size(obj.saccadeIndex,2)-1
                if abs(obj.saccadeIndex(1,a) - obj.saccadeIndex(1,a+1)) >= (floor(obj.stimTime / 5) + 10) % 5ms resolution from DOVES
                    obj.fixationIndex(obj.fixationIndexAdder,1:floor((obj.stimTime / 5))) = obj.saccadeIndex(1,a) + 5 : obj.saccadeIndex(1,a) + (4 + floor(obj.stimTime / 5));
                    obj.fixationIndexAdder = obj.fixationIndexAdder + 1;
                end
            end
            if obj.fixationIndexAdder == 1
                error('No DOVES fixations of suitable length. Please choose shorter stimTime.')
            end
                     
            % Calculate trajectories
            % 1st dimension: different fixational paths
            % 2nd dimension: time
            % 3rd dimension: each amplification
            
            obj.xTraj = zeros(size(obj.fixationIndex,1),size(obj.fixationIndex,2),size(obj.amplification,2));
            obj.yTraj = zeros(size(obj.fixationIndex,1),size(obj.fixationIndex,2),size(obj.amplification,2));
            for a = 1:size(obj.amplification,2)
                for n = 1:size(obj.fixationIndex,1)
                    for m = 1:size(obj.fixationIndex,2)
                        obj.xTraj(n,m,a) = round(obj.fixMovement.x(1,obj.fixationIndex(n,m))*obj.amplification(1,a));
                        obj.yTraj(n,m,a) = round(obj.fixMovement.y(1,obj.fixationIndex(n,m))*obj.amplification(1,a));
                    end
                    obj.xTraj(n,:,a) = -(obj.xTraj(n,:,a) - obj.xTraj(n,1,a));
                    obj.yTraj(n,:,a) = obj.yTraj(n,:,a) - obj.yTraj(n,1,a); 
                end
            end
            
            obj.timeTraj = (0:size(obj.xTraj,2)-1) ./ 200; % DOVES resolution in sec.
                
            % Adjust image
            img = load(strcat('+edu/+washington/+riekelab/+freedland/+images',obj.imageStyle,'.mat'));
            img = double(img.img);
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            obj.imageMatrixFull = zeros(size(img,1),size(img,2),size(obj.zooming,2));
            
            % 1, 2 dimensions: image
            % 3rd dimension: respective zooms
            
            for a = 1:size(obj.zooming,2)
                obj.tempImage = imresize(img,obj.zooming(1,a));
                if obj.zooming(1,a) < 1
                    obj.tempImage2 = repmat(obj.tempImage,ceil(1/obj.zooming(1,a)),ceil(1/obj.zooming(1,a)));
                    obj.imageMatrixFull(:,:,a) = obj.tempImage2(1:size(img,1),1:size(img,2));
                else
                    obj.imageMatrixFull(:,:,a) = obj.tempImage(1:size(img,1),1:size(img,2));
                end
            end

            % Adjust axes and units for monitor
            obj.xTraj = obj.xTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj = obj.yTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');

            % To avoid adaptation, we randomize the order
            obj.selectionIndex = 1;
            
            % All possible combinations of zooms and amplifications
            obj.zoomAmpPerm = [repelem(1:size(obj.zooming,2),1,size(obj.amplification,2));...
                repmat(1:size(obj.amplification,2),1,size(obj.zooming,2))];
            obj.randomIndex = randperm(size(obj.zoomAmpPerm,2));
            
            if obj.randomFixations == true
                obj.randomFixationIndex = randperm(size(obj.xTraj,1));
                obj.fixationIndex = 1;
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.controlTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            % Zoom is the first row
            epoch.addParameter('specificZoom', obj.zooming(obj.zoomAmpPerm(1,obj.randomIndex(obj.selectionIndex))));
            
            % Amplification is the second row
            epoch.addParameter('specificAmplification', obj.amplification(obj.zoomAmpPerm(2,obj.randomIndex(obj.selectionIndex))));
            epoch.addParameter('currentStimSet', obj.currentStimSet);
        end
        
        function p = createPresentation(obj)
            
            % Prep stage for presentation
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.controlTime + obj.tailTime) * 1e-3); % add before/after image as control
            
            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            obj.picture = obj.imageMatrixFull(:,:,obj.zoomAmpPerm(1,obj.randomIndex(obj.selectionIndex)))*255;
            
            if obj.mirroring == true
                obj.mirroredImg = [flip(flip(obj.picture,2)) flip(obj.picture) flip(flip(obj.picture,2));
                    flip(obj.picture,2) obj.picture flip(obj.picture,2); flip(flip(obj.picture,2)) flip(obj.picture) flip(flip(obj.picture,2))];
                obj.imageMatrix = uint8(obj.mirroredImg);
            else
                obj.imageMatrix = uint8(obj.picture);
            end
            
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
            if obj.randomFixations == true
                obj.xTrajAmp = obj.xTraj(obj.randomFixationIndex(obj.fixationIndex),:,obj.zoomAmpPerm(2,obj.randomIndex(obj.selectionIndex)));
                obj.yTrajAmp = obj.yTraj(obj.randomFixationIndex(obj.fixationIndex),:,obj.zoomAmpPerm(2,obj.randomIndex(obj.selectionIndex)));
            else
                obj.xTrajAmp = obj.xTraj(1,:,obj.zoomAmpPerm(2,obj.randomIndex(obj.selectionIndex)));
                obj.yTrajAmp = obj.yTraj(1,:,obj.zoomAmpPerm(2,obj.randomIndex(obj.selectionIndex)));
            end
            
            % Apply eye trajectories to move image around.
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - (obj.preTime+obj.controlTime/2)/1e3, p0));
            
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
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime + obj.controlTime) * 1e-3);
            p.addController(sceneVisible);
            
                % If we are still in one runthru,
            if obj.selectionIndex < size(obj.amplification,2) * size(obj.zooming,2)
                obj.selectionIndex = obj.selectionIndex + 1; % choose next index
                
                % If we reach end of one runthru,
            elseif obj.selectionIndex == size(obj.amplification,2) && obj.numEpochsPrepared < obj.numberOfAverages * size(obj.amplification,2) * size(obj.zooming,2)
                
                % Randomize amplifications/zooms and start again
                obj.selectionIndex = 1;
                obj.randomIndex = randperm(size(obj.amplification,2));
                
                % Randomize fixation path (if desired).
                if obj.randomFixations == true && obj.fixationIndex < floor(obj.numEpochsPrepared / obj.numberOfAverages)
                    obj.fixationIndex = obj.fixationIndex + 1;
                elseif obj.fixationIndex >= floor(obj.numEpochsPrepared / obj.numberOfAverages)
                    obj.randomFixationIndex = randperm(size(obj.xTraj,1));
                    obj.fixationIndex = 1;
                end
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * size(obj.amplification,2) * size(obj.zooming,2);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * size(obj.amplification,2) * size(obj.zooming,2);
        end

    end
    
end