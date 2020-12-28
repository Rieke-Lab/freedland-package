% Reduces the dimensionality of natural images over an array of subunits.
% By J. Freedland, 2020.
classdef reduceNaturalImages < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % Cell information
        rfSigmaCenter = 30; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Subunit information
        subunitLocation_r       = [0 10 20 30];  % grid of subunits (in radial coordinates: r, in microns)
        subunitLocation_theta   = [0 180 0 180]; % grid of subunits (in radial coordinates: theta, in degrees)
        subunitRadius           = 40;            % in microns
        
        % Natural image information
        imageNo = 5;                % natural image to reduce (1 - 101)
        observerNo = 1;             % different observer's eye trajectory (1 - 19)
        showNaturalMovie = true;    % show natural image movie prior to reduced movie

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(3) % number of repeats
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        backgroundIntensity
        movieExport
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
            
            % Load generic settings
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
            
            % Use DOVES eye movements to create movie. Convolve movie with receptive field structure.
            [convolvedStimulus, obj.movieExport{1,1}] = edu.washington.riekelab.freedland.videoGeneration.utils.weightedTrajectory(settings, img, RFFilter);
            
            % The DOVES database is in units of arcmin
            subunitRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            subunitLocationPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitLocation_r,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            
            % Build reduced stimulus
            subunitMasks = zeros([settings.videoSize,length(obj.subunitLocation_r)]);
            reducedStimulus = zeros(size(convolvedStimulus));
            for a = 1:length(obj.subunitLocation_r)
                
                % Build a mask for each subunit
                xCoord = settings.videoSize(2)/2 + subunitLocationPix(a) .* cos(deg2rad(obj.subunitLocation_theta(a)));
                yCoord = settings.videoSize(1)/2 - subunitLocationPix(a) .* sin(deg2rad(obj.subunitLocation_theta(a)));
                [xx,yy] = meshgrid(1:settings.videoSize(2),1:settings.videoSize(1));
                subunitMasks(:,:,a) = sqrt((xx - xCoord).^2 + (yy - yCoord).^2) <= subunitRadiusPix;
                
                % Integrate each frame using linear-equivalent disk
                normalizingValue = sum(RFFilter .* subunitMasks(:,:,a),[1 2]); % RF * mask
                for b = 1:size(obj.stimulus.raw,4)
                    rawValue = sum(convolvedStimulus(:,:,1,b) .* subunitMasks(:,:,a),[1 2]); % RF * image * mask
                    reducedStimulus(:,:,1,b) = reducedStimulus(:,:,1,b) + subunitMasks(:,:,a) .* (rawValue / normalizingValue);
                end
            end
            
            % For overlapping subunits: average results
            obj.movieExport{1,2} = reducedStimulus ./ repmat(sum(subunitMasks,3),1,1,1,size(reducedStimulus,4));
            
            obj.counter = 0;
            if obj.showNaturalMovie == true
                obj.order = repmat([1 2],1,length(obj.numberOfAverages));
            else
                obj.order = repelem(2,length(obj.numberOfAverages));
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            if obj.order(obj.counter+1) == 1
                epoch.addParameter('movieType', 'natural');
            elseif obj.order(obj.counter+1) == 2
                epoch.addParameter('movieType', 'reduced');
            end
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); 
            refreshRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity
            
            % Load individual movie
            individualMovie = obj.movieExport{1,obj.order(obj.counter+1)};
            
            % Prep to display image
            initMatrix = uint8(255 .* obj.backgroundIntensity .* ones(size(obj.stimulus.reduced,1),size(obj.stimulus.reduced,2)));
            mov = stage.builtin.stimuli.Image(initMatrix); % Starting image (blank)
            mov.size = canvasSize;
            mov.position = canvasSize/2;
            mov.setMinFunction(GL.LINEAR); % Can also use GL.NEAREST for no pixel interpolation
            mov.setMagFunction(GL.LINEAR);
            p.addStimulus(mov);
            
            % Program each frame
            preFrames = round(refreshRate * (obj.preTime/1e3));
            movController = stage.builtin.controllers.PropertyController(mov, 'imageMatrix',...
                @(state)getMovieFrame(obj, state.frame - preFrames, individualMovie));
            p.addController(movController);
            
            function frame = getMovieFrame(obj, frame, movie)
                if frame < 0 % Before movie begins
                    frame = uint8(obj.backgroundIntensity .* 255 .* ones(size(movie,1),size(movie,2)));
                else
                    frame = uint8(movie(:,:,1,frame));
                end
            end

            % Hide during pre & post
            movVisible = stage.builtin.controllers.PropertyController(mov, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(movVisible); 
            
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.order);
        end
    end
end