% Reduces the dimensionality of natural images over an array of subunits.
% By J. Freedland, 2020.
classdef reduceNaturalImages < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime     = 250  % in ms
        stimTime    = 5500 % in ms
        tailTime    = 250  % in ms
        
        % Cell information
        rfSigmaCenter   = 30;  % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Subunit information
        subunitDiameter     = 40;          % in microns
        subunitLocation_y   = [25 50 75];  % grid of subunits (in um, first row from exported subunits.txt)
        subunitLocation_x   = [0 120 240]; % grid of subunits (in um, second row from exported subunits.txt)
        
        % Natural image information
        imageNo     = 5; % natural image to reduce (1 - 101)
        observerNo  = 1; % different observer's eye trajectory (1 - 19)
        
        % Experiment settings
        experimentSettings       = 'add subunits';  % "simple": show specified subunits; "add subunits": show groups of (1, 2, 3, ...) subunits.
        showRandomCoordinates    = true;            % include trials with randomly-located subunits. Follows same approach as "experimentSettings".
        showNaturalMovie         = true;            % show natural image movie prior to reduced movie
        naturalMovieRegion       = 'center-only';   % specific region to show natural image movie
        showLinearEquivalentDisk = true;            % show reduced movie (single linear equivalent disk)
        randomize                = true;            % display selected movies in random order

        % Additional parameters
        onlineAnalysis      = 'extracellular'
        numberOfAverages    = uint16(3) % number of repeats
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        experimentSettingsType = symphonyui.core.PropertyType('char', 'row', {'simple', 'add subunits'}) 
        naturalMovieRegionType = symphonyui.core.PropertyType('char', 'row', {'center-only', 'surround-only', 'full-field'}) 
        backgroundIntensity
        randomCoordinates_x
        randomCoordinates_y
        directory
        filename
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
            [convolvedStimulus, rawMov] = edu.washington.riekelab.freedland.videoGeneration.utils.weightedTrajectory(settings, img, RFFilter);
            
            % Adjust units to DOVES database (units of arcmin)
            subunitRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitDiameter/2,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            subunitLocationPix_x = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitLocation_x,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            subunitLocationPix_y = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitLocation_y,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            
            % Identify natural image movie region to show.
            [xx,yy] = meshgrid(1:settings.videoSize(2),1:settings.videoSize(1));
            m = sqrt((xx - settings.videoSize(2)/2).^2 + (yy - settings.videoSize(1)/2).^2);

            % Refine naturalistic movie
            if obj.showNaturalMovie == true
                if strcmp(obj.naturalMovieRegion,'center-only')
                    mask = m <= settings.rfSizing.zeroPt;
                elseif strcmp(obj.naturalMovieRegion,'surround-only')
                    mask = m >= settings.rfSizing.zeroPt & m < max(settings.diskRadii);
                elseif strcmp(obj.naturalMovieRegion,'full-field')
                    mask = m < max(settings.diskRadii);
                end
                
                % Only display for desired spatial region.
                for a = 1:size(rawMov,4)
                    rawMov(:,:,1,a) = rawMov(:,:,1,a) .* mask + ~mask .* obj.backgroundIntensity;
                end
                obj.filename{1,1} = 'natural';
                outputMovie{1,1} = rawMov .* 255;
            end
            
            % Calculate linear equivalent center
            if obj.showLinearEquivalentDisk == true
                linearEquiv = zeros(size(rawMov));
                mask = m <= settings.rfSizing.zeroPt; % Identify center
                for a = 1:size(convolvedStimulus,4)
                    rawValue = convolvedStimulus(:,:,1,a) .* mask;
                    rawValue = sum(rawValue(:));
                    
                    normalizingValue = RFFilter .* mask;
                    normalizingValue = sum(normalizingValue(:));
                    
                    linearEquiv(:,:,1,a) = mask .* (rawValue ./ normalizingValue) + ~mask .* obj.backgroundIntensity;
                end
                obj.filename{2,1} = 'linearEquivalent';
                outputMovie{2,1} = linearEquiv .* 255;
            end
            
            %%% Build reduced stimulus
            subunitMasks = zeros([settings.videoSize,length(subunitLocationPix_x)]);
            [xx,yy] = meshgrid(1:settings.videoSize(2),1:settings.videoSize(1));
            
            % Generate random coordinates as needed
            if obj.showRandomCoordinates == true
                obj.randomCoordinates_x = (rand(length(subunitLocationPix_x),1) - 0.5) .* settings.rfSizing.zeroPt;
                obj.randomCoordinates_y = (rand(length(subunitLocationPix_y),1) - 0.5) .* settings.rfSizing.zeroPt;
                
                % Confirm coordinates are sufficiently far apart
                for a = 1:length(obj.randomCoordinates_x)-1
                    referenceCoordinate = [obj.randomCoordinates_x(a), obj.randomCoordinates_y(a)];
                    for b = a+1:length(obj.randomCoordinates_x)
                        iter = 1;
                        while (norm(referenceCoordinate - [obj.randomCoordinates_x(b), obj.randomCoordinates_y(b)]) < subunitRadiusPix*2) && (iter < 1e4)
                            obj.randomCoordinates_x = (rand(length(subunitLocationPix_x),1) - 0.5) .* settings.rfSizing.zeroPt;
                            obj.randomCoordinates_y = (rand(length(subunitLocationPix_x),1) - 0.5) .* settings.rfSizing.zeroPt;
                            iter = iter+1;
                        end
                    end
                end
                iterations = 2;
            else
                iterations = 1;
            end
            
            for d = 1:iterations % Experiment set
                
                % Whether to produce one movie or multiple movies
                if strcmp(obj.experimentSettings,'simple')
                    reducedStimulus = zeros(size(convolvedStimulus));
                elseif strcmp(obj.experimentSettings,'add subunits')
                    reducedStimulus = zeros([size(convolvedStimulus),length(subunitLocationPix_x)]);
                end
                
                for a = 1:length(subunitLocationPix_x)
                    
                    % Build a mask for each subunit
                    if d == 1 % Desired coordinates
                        xCoord = settings.videoSize(2)/2 + subunitLocationPix_y(a); % intentionally flipped (meshgrid rotates)
                        yCoord = settings.videoSize(1)/2 + subunitLocationPix_x(a);
                    elseif d == 2 % Random coordinates
                        xCoord = settings.videoSize(2)/2 + obj.randomCoordinates_y(a); % intentionally flipped (meshgrid rotates)
                        yCoord = settings.videoSize(1)/2 + obj.randomCoordinates_x(a);
                    end
                    subunitMask = sqrt((xx - xCoord).^2 + (yy - yCoord).^2) <= subunitRadiusPix;

                    % Integrate each frame using linear-equivalent disk
                    normalizingValue = RFFilter .* subunitMask;
                    normalizingValue = sum(normalizingValue(:)); % RF * mask
                    for b = 1:size(rawMov,4)
                        rawValue = convolvedStimulus(:,:,1,b) .* subunitMask;
                        rawValue = sum(rawValue(:)); % RF * image * mask

                        % Add each frame
                        if strcmp(obj.experimentSettings,'simple')
                            reducedStimulus(:,:,1,b) = reducedStimulus(:,:,1,b) + subunitMask .* (rawValue / normalizingValue);
                        elseif strcmp(obj.experimentSettings,'add subunits')
                            for c = a:length(subunitLocationPix_x)
                                reducedStimulus(:,:,1,b,c) = reducedStimulus(:,:,1,b,c) + subunitMask .* (rawValue / normalizingValue);
                            end
                        end
                    end
                end
                reducedMov = reducedStimulus ./ repmat(sum(subunitMasks,3),1,1,1,size(reducedStimulus,4),size(reducedStimulus,5)); % For overlapping subunits: average results
                reducedMov(isnan(reducedMov)) = obj.backgroundIntensity;

                for a = 1:size(reducedMov,5)
                    if d == 1
                        obj.filename = [obj.filename;{strcat('reduced_',mat2str(a),'_',mat2str([subunitLocationPix_x(a) subunitLocationPix_y(a)]))}];
                    elseif d == 2
                        obj.filename = [obj.filename;{strcat('reducedRandomized_',mat2str(a),'_',mat2str([obj.randomCoordinates_x(a) obj.randomCoordinates_y(a)]))}];
                    end
                    outputMovie = [outputMovie;{reducedMov .* 255}];
                end
            end
            
            % Export all movies as .avi (no compression)
            obj.directory = 'Documents/freedland-package/+edu/+washington/+riekelab/+freedland/+movies/';
            
            % Add pre-time and post-time manually
            refreshRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
            blankFrames = ones(size(reducedMov(:,:,1,1))) .* 255 .* obj.backgroundIntensity;
            preFrames = repmat(blankFrames,1,1,1,round(refreshRate * (obj.preTime/1e3)));
            postFrames = repmat(blankFrames,1,1,1,round(refreshRate * (obj.tailTime/1e3)));

            for a = 1:length(outputMovie)
                if ~isempty(outputMovie{a,1})
                    specificMov = uint8(cat(4,preFrames,outputMovie{a,1},postFrames));
                    
                    % Export movies
                    v = VideoWriter(strcat(obj.directory,obj.filename{a,1}),'Uncompressed AVI');
                    v.FrameRate = refreshRate;
                    
                    open(v)
                    for b = 1:size(specificMov,4)
                        writeVideo(v,specificMov(:,:,1,b));
                    end
                    close(v)
                end
            end
            
            obj.counter = 0;
            obj.order = repmat(obj.filename,obj.numberOfAverages,1);
            if obj.randomize == true
                obj.order = obj.order(randperm(length(obj.order)));
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('movieType', obj.filename{obj.order(obj.counter+1)});
            
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
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity

            % Prep to display image
            scene = stage.builtin.stimuli.Movie(fullfile(obj.directory,strcat(obj.order{obj.counter+1},'.avi')));
            scene.size = canvasSize;
            p0 = canvasSize/2;
            scene.position = p0;
            scene.setMinFunction(GL.LINEAR); % Linear scaling to monitor
            scene.setMagFunction(GL.LINEAR);
            p.addStimulus(scene);
            
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