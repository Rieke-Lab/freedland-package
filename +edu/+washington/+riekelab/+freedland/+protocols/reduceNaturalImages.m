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
        subunitLocation_x   = [10 40 -30];  % grid of subunits (in um, first row from exported subunits.txt)
        subunitLocation_y   = [0 -30 -40]; % grid of subunits (in um, second row from exported subunits.txt)
        
        % Natural image information
        imageNo     = 5; % natural image to reduce (1 - 101)
        observerNo  = 1; % different observer's eye trajectory (1 - 19)
        
        % Experiment settings
        experimentSettings       = 'add subunits';  % "simple": show specified subunits; "add subunits": show groups of (1, 2, 3, ...) subunits.
        showRandomCoordinates    = true;            % include trials with randomly-located subunits. Follows same approach as "experimentSettings".
        showNaturalMovie         = true;            % show natural image movie prior to reduced movie
        naturalMovieRegion       = 'center-only';   % specific region to show natural image movie
        showLinearEquivalentDisk = true;            % show reduced movie (single linear equivalent disk)
        showNaturalDots          = true;            % show sets of movies with only dots containing natural movie.
        randomize                = true;            % display selected movies in random order

        % Additional parameters
        onlineAnalysis      = 'extracellular'
        numberOfAverages    = uint16(5) % number of repeats
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        experimentSettingsType = symphonyui.core.PropertyType('char', 'row', {'simple', 'add subunits'}) 
        naturalMovieRegionType = symphonyui.core.PropertyType('char', 'row', {'center-only', 'surround-only', 'full-field'}) 
        backgroundIntensity
        coordinates
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

            outputMovie = {};
            obj.filename = {};
            obj.coordinates = {};
            
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
                obj.coordinates{1,1} = obj.naturalMovieRegion; % x
                obj.coordinates{1,2} = obj.naturalMovieRegion; % y
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
                obj.coordinates{2,1} = 'linearEquivalent'; % x
                obj.coordinates{2,2} = 'linearEquivalent'; % y
            end
            
            %%% Build reduced stimulus
            [xx,yy] = meshgrid(1:settings.videoSize(2),1:settings.videoSize(1));
            
            % Generate random coordinates as needed (in arcmin)
            if obj.showRandomCoordinates == true
                randomCoordinates_x = (zeros(length(subunitLocationPix_x),1));
                randomCoordinates_y = (zeros(length(subunitLocationPix_x),1));
                
                for a = 1:length(randomCoordinates_x)
                    % Dot cannot exceed center
                    randomCoordinates_x(a) = (rand() - 0.5) .* 2 .* (settings.rfSizing.zeroPt - subunitRadiusPix);
                    randomCoordinates_y(a) = (rand() - 0.5) .* 2 .* (settings.rfSizing.zeroPt - subunitRadiusPix);
                    
                    % Confirm subunit is sufficiently far away from others
                    A = repmat([randomCoordinates_x(a) randomCoordinates_y(a)],a-1,1);
                    B = [randomCoordinates_x(1:a-1) randomCoordinates_y(1:a-1)];
                    if ~isempty(B)
                        iterChecker = 1;
                        while sum(sqrt(sum((A - B).^2,2)) < subunitRadiusPix*2) > 0 % Too close
                            randomCoordinates_x(a) = (rand() - 0.5) .* 2 .* (settings.rfSizing.zeroPt - subunitRadiusPix);
                            randomCoordinates_y(a) = (rand() - 0.5) .* 2 .* (settings.rfSizing.zeroPt - subunitRadiusPix);
                            A = repmat([randomCoordinates_x(a) randomCoordinates_y(a)],a-1,1);
                            B = [randomCoordinates_x(1:a-1) randomCoordinates_y(1:a-1)];
                            
                            iterChecker = iterChecker + 1;
                            if iterChecker > 1e4
                                error('Unable to find reasonable random coordinates.')
                            end
                        end
                    end
                end

                % For fairness: order by most-likely subunits (center)
                [~,i] = sort(randomCoordinates_x.^2 + randomCoordinates_y.^2);
                randomCoordinates_x = randomCoordinates_x(i);
                randomCoordinates_y = randomCoordinates_y(i);
                iterations = 2;
            else
                iterations = 1;
            end
            
            tic
            % Build reduced representations
            for d = 1:iterations
                subunitMaskTracker = zeros([size(xx) length(subunitLocationPix_x)]);
                
                % Whether to produce one movie or multiple movies
                if strcmp(obj.experimentSettings,'simple')
                    reducedStimulus = zeros(size(convolvedStimulus));
                elseif strcmp(obj.experimentSettings,'add subunits')
                    reducedStimulus = zeros([size(convolvedStimulus),length(subunitLocationPix_x)]);
                end
                
                if obj.showNaturalDots == true
                    naturalReducedStimulus = reducedStimulus;
                end
                
                for a = 1:length(subunitLocationPix_x)
                    % Build a mask for each subunit
                    if d == 1 % Desired coordinates
                        xCoord = settings.videoSize(2)/2 + subunitLocationPix_x(a); % intentionally flipped (meshgrid rotates)
                        yCoord = settings.videoSize(1)/2 + subunitLocationPix_y(a);
                    elseif d == 2 % Random coordinates
                        xCoord = settings.videoSize(2)/2 + randomCoordinates_x(a); % intentionally flipped (meshgrid rotates)
                        yCoord = settings.videoSize(1)/2 + randomCoordinates_y(a);
                    end
                    subunitMask = sqrt((xx - xCoord).^2 + (yy - yCoord).^2) <= subunitRadiusPix;
                    subunitMaskTracker(:,:,a) = subunitMask;

                    % Integrate each frame using linear-equivalent disk
                    normalizingValue = RFFilter .* subunitMask;
                    normalizingValue = sum(normalizingValue(:)); % RF * mask
                    for b = 1:size(rawMov,4)
                        rawValue = convolvedStimulus(:,:,1,b) .* subunitMask;
                        rawValue = sum(rawValue(:)); % RF * image * mask

                        % Add each frame
                        if strcmp(obj.experimentSettings,'simple')
                            reducedStimulus(:,:,1,b) = reducedStimulus(:,:,1,b) + subunitMask .* (rawValue / normalizingValue);
                            if obj.showNaturalDots == true
                                naturalReducedStimulus(:,:,1,b) = naturalReducedStimulus(:,:,1,b) + subunitMask .* rawMov(:,:,1,b);
                            end
                        elseif strcmp(obj.experimentSettings,'add subunits')
                            for c = a:length(subunitLocationPix_x)
                                reducedStimulus(:,:,1,b,c) = reducedStimulus(:,:,1,b,c) + subunitMask .* (rawValue / normalizingValue);
                                if obj.showNaturalDots == true
                                    naturalReducedStimulus(:,:,1,b,c) = naturalReducedStimulus(:,:,1,b,c) + subunitMask .* rawMov(:,:,1,b);
                                end
                            end
                        end
                    end
                end
                
                % Average overlapping regions
                for c = 1:length(subunitLocationPix_x)
                    reducedStimulus(:,:,:,:,c) = reducedStimulus(:,:,:,:,c) ./ repmat(sum(subunitMaskTracker(:,:,1:c),3),1,1,1,size(reducedStimulus,4)); % For overlapping subunits: average results
                    if obj.showNaturalDots == true
                        naturalReducedStimulus(:,:,:,:,c) = naturalReducedStimulus(:,:,:,:,c) ./ repmat(sum(subunitMaskTracker(:,:,1:c),3),1,1,1,size(naturalReducedStimulus,4));
                    end
                end
                reducedStimulus(isnan(reducedStimulus)) = obj.backgroundIntensity;
                if obj.showNaturalDots == true
                    naturalReducedStimulus(isnan(naturalReducedStimulus)) = obj.backgroundIntensity;
                end

                for a = 1:size(reducedStimulus,5)
                    if d == 1
                        if strcmp(obj.experimentSettings,'simple')
                            C = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                            [subunitLocationPix_x; subunitLocationPix_y],obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2um'),1);
                            specificFilename = 'reduced_all';
                        elseif strcmp(obj.experimentSettings,'add subunits') 
                            C = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                            [subunitLocationPix_x(1:a); subunitLocationPix_y(1:a)],obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2um'),1);
                            specificFilename = strcat('reduced_',mat2str(a),'-subunits'); % Coordinates in microns
                        end
                    elseif d == 2
                        if strcmp(obj.experimentSettings,'simple')
                            C = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                            [randomCoordinates_x'; randomCoordinates_y'],obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2um'),1);
                            specificFilename = 'reducedRandomized_all'; % Coordinates in microns
                        elseif strcmp(obj.experimentSettings,'add subunits')    
                            C = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                            [randomCoordinates_x(1:a)'; randomCoordinates_y(1:a)'],obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2um'),1);
                            specificFilename = strcat('reducedRandomized_',mat2str(a),'-subunits'); % Coordinates in microns
                        end
                    end
                    obj.filename = [obj.filename;{specificFilename}]; % Coordinates in microns
                    obj.coordinates = cat(1,obj.coordinates,[{C(1,:)},{C(2,:)}]); %[ x y ]
                    outputMovie = [outputMovie;{reducedStimulus(:,:,:,:,a) .* 255}];
                    
                    if obj.showNaturalDots == true
                        obj.coordinates = cat(1,obj.coordinates,[{C(1,:)},{C(2,:)}]); %[ x y ]
                        obj.filename = [obj.filename;{strcat('naturalDot-',specificFilename)}];
                        outputMovie = [outputMovie;{naturalReducedStimulus(:,:,:,:,a) .* 255}];
                    end
                end
            end
            toc
            
            % Export all movies
            obj.directory = 'Documents/freedland-package/+edu/+washington/+riekelab/+freedland/+movies/';
            
            % Add pre-time and post-time manually
            refreshRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
            blankFrames = ones(size(reducedStimulus(:,:,1,1))) .* 255 .* obj.backgroundIntensity;
            preFrames = repmat(blankFrames,1,1,1,round(refreshRate * (obj.preTime/1e3)));
            postFrames = repmat(blankFrames,1,1,1,round(refreshRate * (obj.tailTime/1e3)));

            disp('Exporting movies...')
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
            obj.order = repmat(1:length(obj.filename),1,obj.numberOfAverages);
            if obj.randomize == true
                obj.order = obj.order(randperm(length(obj.order)));
            end
            
            t = length(obj.order) .* ((obj.preTime+obj.stimTime+obj.tailTime)/1000 + 1.5) / 60; % in minutes
            disp(strcat('Estimated stimulus time (+1.5 sec delay): ',mat2str(round(t,1)),'min.'))
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('movieName', obj.filename{obj.order(obj.counter+1)});
            epoch.addParameter('xCoordinate_um', obj.coordinates{obj.order(obj.counter+1),1});
            epoch.addParameter('yCoordinate_um', obj.coordinates{obj.order(obj.counter+1),2});
            
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
            scene = stage.builtin.stimuli.Movie(fullfile(obj.directory,strcat(obj.filename{obj.order(obj.counter+1)},'.avi')));
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