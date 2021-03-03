% Plays a natural image movie, courtesy of the DOVES database.
classdef sparseNaturalImageMovie < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        
        % Stimulus timing
        preTime     = 250   % in ms
        stimTime    = 5500  % in ms
        tailTime    = 250   % in ms
        
        % Cell details
        rfSigmaCenter   = 30;  % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Natural image trajectory
        imageNo     = 5;    % natural image number (1 to 101)
        observerNo  = 1;    % observer number (1 to 19)
        
        % Mask details
        coverageStyle   = 'static';  % cover the center with a mask of a given type.
        areaCoverage    = 60;        % percent area of natural image to hide.
        uniqueCoverage  = 5;         % number of distinct ways to cover natural image.
        randomize       = true;      % whether to randomize presentation of different coverings.
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
        
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        coverageStyleType = symphonyui.core.PropertyType('char', 'row', {'static', 'equivalentDisk'}) 
        backgroundIntensity
        xTraj
        yTraj
        timeTraj
        luminanceTraj
        imageMatrix
        surroundMask
        centerMask
        noiseSeed
        counter
        order
        diskLuminance
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
            
            %%% Natural image movie
            settings = edu.washington.riekelab.freedland.videoGeneration.demoUtils.loadSettings(...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'),...
                obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'),...
                obj.rfSigmaCenter,obj.rfSigmaSurround);
            
            % Pull base trajectories and image information.
            [path,img] = edu.washington.riekelab.freedland.videoGeneration.utils.pathDOVES(obj.imageNo, obj.observerNo);
            
            % Identify image properties
            img = img./max(img(:));
            obj.backgroundIntensity = mean(img(:)); % Background light intensity
            settings.imageMatrix = uint8(img*255);

            % Eye movement patterns from DOVES database.
            DOVES_trajectory    = 1/200:1/200:(length(path.x)/200); % 200 Hz in DOVES database
            monitorTrajectory   = 1/settings.monitorFrameRate:1/settings.monitorFrameRate:(length(path.x)/200);
            settings.xTraj = interp1(DOVES_trajectory,path.x,monitorTrajectory);
            settings.yTraj = interp1(DOVES_trajectory,path.y,monitorTrajectory);

            % Calculate receptive field structure
            [RFFilter,settings.rfSizing] = edu.washington.riekelab.freedland.videoGeneration.rfUtils.calculateFilter(settings);
            
            % Use DOVES eye movements to create movie.
            [~, naturalMovie] = edu.washington.riekelab.freedland.videoGeneration.utils.weightedTrajectory(settings, img, RFFilter); 
            
            % Convert from arcmin to monitor pix
            naturalMovie = imresize(naturalMovie,settings.monitorSize);
            RFFilter = imresize(RFFilter,settings.monitorSize);
            %%%
            
            %%% Retrofit movie for animated monitor
            % Invert for monitor and adjust position relative to center of image
            obj.xTraj = -(path.x - size(img,2)/2);
            obj.yTraj = (path.y - size(img,1)/2);
            obj.imageMatrix = settings.imageMatrix;
            
            % Convert from DOVES units (1 px = 1 arcmin) to monitor units
            obj.xTraj = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.xTraj,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
            obj.yTraj = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.yTraj,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % convert DOVES resolution (200Hz) to seconds
            %%%
            
            % Surround mask
            [xx,yy] = meshgrid(1:settings.monitorSize(2),1:settings.monitorSize(1));
            r = sqrt((xx - settings.monitorSize(2)/2).^2 + (yy - settings.monitorSize(1)/2).^2); 
            centerRadius = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                settings.rfSizing.zeroPt,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
            obj.surroundMask = r > centerRadius;
            
            % Design center masks
            obj.noiseSeed = cell(obj.uniqueCoverage,1);
            obj.centerMask = zeros([size(xx) obj.uniqueCoverage]);
            obj.diskLuminance = zeros(size(naturalMovie,4),obj.uniqueCoverage);
            obj.luminanceTraj = (0:(size(naturalMovie,4)-1)) ./ 60; % convert DOVES resolution (200Hz) to seconds

            for a = 1:obj.uniqueCoverage
                % Generate random seed
                obj.noiseSeed{a} = RandStream.shuffleSeed;
                noiseStream = RandStream('mt19937ar', 'Seed', obj.noiseSeed{a});
                A = noiseStream.rand(size(xx)); % Uniformly distributed [0, 1]
                obj.centerMask(:,:,a) = A < (obj.areaCoverage/100) .* ~obj.surroundMask;

                % Calculate linear equivalency
                if strcmp(obj.coverageStyle,'equivalentDisk')
                    normValue = obj.centerMask(:,:,a) .* RFFilter;
                    normValue = sum(normValue(:));
                    for b = 1:size(naturalMovie,4)
                        rawValue = naturalMovie(:,:,:,b) .* obj.centerMask(:,:,a) .* RFFilter;
                        rawValue = sum(rawValue(:));

                        obj.diskLuminance(b,a) = rawValue / normValue;
                    end
                end
            end
            
            % Add null case for linear equivalent disk
            obj.noiseSeed = [{'natural'};obj.noiseSeed];
            obj.centerMask = cat(3,zeros(size(xx)),obj.centerMask);
            obj.diskLuminance = [zeros(size(naturalMovie,4),1) obj.diskLuminance];
            
            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(length(obj.noiseSeed));
            else
                obj.order = 1:length(obj.noiseSeed);
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;

            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('noiseSeed', obj.noiseSeed{obj.order(obj.counter+1)});
            
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
            
            % Add information to Stage
            p.addStimulus(scene);
            p.addController(scenePosition);
            
            % Add additional controller: mediate when stimulus is visible.
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            % Apply surround mask
            surround = stage.builtin.stimuli.Rectangle();
            surround.position = canvasSize/2;
            surround.size = [canvasSize(1) canvasSize(2)];
            surround.color = obj.backgroundIntensity;
            annulusS = uint8(obj.surroundMask*255);
            surroundMaskX = stage.core.Mask(annulusS);
            surround.setMask(surroundMaskX);
            p.addStimulus(surround);
            
            if ~strcmp(obj.noiseSeed{obj.order(obj.counter+1)},'natural')
                specificCenter = obj.centerMask(:,:,obj.order(obj.counter+1));
                specificLuminance = obj.diskLuminance(:,obj.order(obj.counter+1));

                center = stage.builtin.stimuli.Rectangle();
                center.position = canvasSize/2;
                center.size = [canvasSize(1) canvasSize(2)];
                centerS = uint8(specificCenter*255);
                centerMaskX = stage.core.Mask(centerS);
                center.setMask(centerMaskX);

                if strcmp(obj.coverageStyle,'static')
                    center.color = obj.backgroundIntensity;
                    p.addStimulus(center); % Add stimulus
                elseif strcmp(obj.coverageStyle,'equivalentDisk')
                    centerLED = stage.builtin.controllers.PropertyController(center,...
                    'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, specificLuminance));
                    
                    centerVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                    
                    p.addStimulus(center); % Add stimulus
                    p.addController(centerVisible);
                    p.addController(centerLED)
                end
            end
            
            obj.counter = mod(obj.counter + 1,length(obj.order));
            
            % Match mask color to linear equivalency in real time.
            function s = getBackground(obj, time, equivalency)
                if time < 0
                    s = obj.backgroundIntensity;
                elseif time > obj.timeTraj(end)
                    s = equivalency(end); % Hold last color
                else
                    s = interp1(obj.luminanceTraj,equivalency,time);
                end
            end
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.order) .* obj.numberOfAverages;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.order) .* obj.numberOfAverages;
        end
    end
end