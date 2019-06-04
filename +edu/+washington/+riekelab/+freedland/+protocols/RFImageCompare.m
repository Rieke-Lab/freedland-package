% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFImageCompare < edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        imageNo = 1; % natural image number (1 to 101)
        observerNo = 1; % observer number (1 to 19)
        amplification = 1; % amplify fixations by X. Setting to 0 produces a saccade-only trajectory. 
        trajectory = 'both'; % which type of stimulus to present: natural image, filtered image, or both?        
        filteredImage = 'raw'; % type of filtered image to present
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(10) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        trajectoryType = symphonyui.core.PropertyType('char', 'row', {'natural','filtered','both'})
        filteredImageType = symphonyui.core.PropertyType('char', 'row', {'raw','complement','negative complement'})
        backgroundIntensity
        imageMatrix
        imageMatrix2
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)

            prepareRun@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj);

            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            if strcmp(obj.trajectory,'both') % Splits the epoch into two for online analysis.
                obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',2);
            else % Keeps as a single epoch.
                obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            end
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Pull base trajectories and image information.
            [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.imageNo, obj.observerNo,...
                    'amplification', obj.amplification,'mirroring', true);
                
            % Scale pixels in image to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix = uint8(img2);
            
            % load alternative image
            imageString = num2str(obj.imageNo);
            directory = '+edu/+washington/+riekelab/+freedland/+images/';

            while size(imageString,2) < 3
                imageString = strcat('0',imageString);
            end
            if ip.Results.directory
                filename = strcat(directory,'altimg',imageString,'.mat');
            else
                filename = strcat('altimg',imageString,'.mat');
            end

            load(filename)
            
            if strcmp(obj.filteredImage,'raw')
                B = uint8(picture); % pre-filtered
            elseif strcmp(obj.filteredImage,'complement')
                B = img2 - picture; % relative to original image
            else
                B = picture - img2; % relative to filtered image
            end
            
            obj.imageMatrix2 = uint8(B);
            
            % Produce trajectories
            obj.xTraj = baseMovement.x + fixMovement.x;
            obj.yTraj = baseMovement.y + fixMovement.y;
            
            % Adjust axes and units for monitor (VH Pixels)
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            obj.xTraj = obj.xTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj = obj.yTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            if ~strcmp(obj.trajectory,'both')
                duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            else
                duration = 2 * (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            end
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
                                    
            if ~strcmp(obj.trajectory,'both')
                p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            else % Present both stimuli in succession
                p = stage.core.Presentation(2 * (obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            end

            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(obj.imageMatrix);
            scene2 = stage.builtin.stimuli.Image(obj.imageMatrix2);
            
            scene.size = [size(obj.imageMatrix,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.imageMatrix,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            scene2.size = scene.size;
            p0 = canvasSize/2;
            scene.position = p0;
            scene2.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            
            scene2.setMinFunction(GL.LINEAR);
            scene2.setMagFunction(GL.LINEAR);
            
            % Apply eye trajectories to move image around.
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3, p0));
            scene2Position = stage.builtin.controllers.PropertyController(scene2,...
                'position', @(state)getScenePosition(obj, state.time - cycleType - obj.preTime/1e3, p0));
            
            function p = getScenePosition(obj, time, p0)
                if time <= 0
                    p = p0;
                elseif time > obj.timeTraj(end) % Beyond eye trajectory, hang on last frame
                    p(1) = p0(1) + obj.xTraj(end);
                    p(2) = p0(2) + obj.yTraj(end);
                else % Within eye trajectory and stim time
                    dx = interp1(obj.timeTraj,obj.xTraj,time);
                    dy = interp1(obj.timeTraj,obj.yTraj,time);
                    p(1) = p0(1) + dx;
                    p(2) = p0(2) + dy;
                end
            end
            
            p.addStimulus(scene);
            p.addStimulus(scene2);
            p.addController(scenePosition);
            p.addController(scene2Position);

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3); % 1st
            scene2Visible = stage.builtin.controllers.PropertyController(scene, 'visible', ... % 2nd
                @(state)state.time >= cycleTime + obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime + cycleTime) * 1e-3);
            p.addController(sceneVisible);
            p.addController(scene2Visible);
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages;
        end
    end
end