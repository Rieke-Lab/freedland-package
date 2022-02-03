% Plays a natural image movie, courtesy of the DOVES database.
classdef blurNaturalImageMovie < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        
        % Stimulus timing
        preTime     = 250   % in ms
        stimTime    = 5500  % in ms
        tailTime    = 250   % in ms
        
        % Natural image trajectory
        imageNo     = 5;    % natural image number (1 to 101)
        observerNo  = 1;    % observer number (1 to 19)
        
        % Mask information
        blurSigma     = 200; % sigma of Gaussian blur (in microns)
        maskDiameter  = 300; % only present natural image in RF center (in microns). Set to 0 to ignore.

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(1) % number of epochs to queue
        amp % Output amplifier
        
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        backgroundIntensity
        xTraj
        yTraj
        timeTraj
        imageMatrix
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
            
            % Gather natural image information
            [path,image] = edu.washington.riekelab.freedland.videoGeneration.utils.pathDOVES(obj.imageNo, obj.observerNo);
            image = image ./ max(image(:)); % Normalize (scale to monitor)
            obj.backgroundIntensity = mean(image(:));
            obj.imageMatrix = uint8(image.*255);
                    
            % Isolate DOVES eye trajectories
            obj.xTraj = path.x;
            obj.yTraj = path.y;
            
            % Invert for monitor and adjust position relative to center of image
            obj.xTraj = -(obj.xTraj - size(image,2)/2);
            obj.yTraj = (obj.yTraj - size(image,1)/2);
            
            % Convert from DOVES units (1 px = 1 arcmin) to monitor units
            obj.xTraj = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.xTraj,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
            obj.yTraj = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.yTraj,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % convert DOVES resolution (200Hz) to seconds
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;

            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            % Add metadata from Stage.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % in normal pixels            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity)
            
            % Insert image and sizing information for stage.
            scene = stage.builtin.stimuli.Image(obj.imageMatrix);
            scene.size = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(fliplr(size(obj.imageMatrix)),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix'); % Convert to monitor units
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            % Control the position of the image as a function of time.
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3, p0));

            function p = getScenePosition(obj, time, p0)
                if time <= 0 % Before stimulus begins
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
            
            %%% Create Gaussian envelope (pasted from Stage's code)
            resolution = 512; % Stage default
            step = 2 / (resolution - 1);
            [xx, yy] = meshgrid(-1:step:1, -1:step:1);
            distanceMatrix = sqrt(xx.^2 + yy.^2);

            % Convert from um to canvas pixels
            sigma = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(obj.blurSigma,...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            gaussian = uint8(exp(-distanceMatrix.^2 / (2 * sigma^2)) * 255);
            scene.setMask = stage.core.Mask(gaussian);
            %%%
            
            % Add information to Stage
            p.addStimulus(scene);
            p.addController(scenePosition);
            
            % Add additional controller: mediate when stimulus is visible.
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            if (obj.maskDiameter > 0) % Create mask
                
                maskDiameterPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(obj.maskDiameter,...
                    obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
                
                mask = stage.builtin.stimuli.Ellipse();
                mask.position = canvasSize/2;
                mask.color = obj.backgroundIntensity;
                mask.radiusX = maskDiameterPix/2;
                mask.radiusY = maskDiameterPix/2;
                p.addStimulus(mask); %add mask
            end
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages;
        end
    end
end
