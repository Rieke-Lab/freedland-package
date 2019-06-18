% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFDiskArrayProjection < edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        imageNo = 5; % natural image number (1 to 101)
        observerNo = 1; % observer number (1 to 19)
        amplification = 1; % amplify fixations by X. Setting to 0 produces a saccade-only trajectory. 
        trajectory = 'both'; % which type of stimulus to present: natural image, disk replacement, or both?        
        
        % RF field information
        rfSigmaCenter = 70; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 170; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Disk placement
        disks = 3; % number of disks that are evenly placed over a natural image.
        overrideRadii = [0 0.75 2 3]; % only takes effect if any value is >0. Allows any number of disks in any distribution. In pixels, for a 800px monitor, must contain 0 and 400 and can be represented as: [0 50 100 400]. In RF coordinates, must contain 0 and 3, where 1, 2 are the radius of the center, surround respectively. For a [center surround] = [70 170], [0 0.5 1 1.5 2 3] in RF = [0 35 70 120 170 400] in px.
        overrideCoordinate = 'RF'; % type of coordinates to measure disk radii.
        xSliceFrequency = 1; % how many radial slices to cut between 0 and 90 degrees.
        ySliceFrequency = 1; % how many radial slices to cut between 90 and 180 degrees.
        rotateSlices = 0; % degrees to rotate all slices (0 to 90).
        disksIgnoreCut = [3 0]; % starting from the center disk and moving outwards, how many disks should we NOT cut (keep circular)?

        % Disk type
        meanDisks = [1 2 3]; % starting from the center disk and moving outwards, which disks should be averaged?
        naturalDisks = [0 0];  % starting from the center disk and moving outwards, which disks should remain a natural image?
        backgroundDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be left at background intensity?
        meanIntegration = 'gaussian';
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(10) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        trajectoryType = symphonyui.core.PropertyType('char', 'row', {'natural','disk','both'})
        meanIntegrationType = symphonyui.core.PropertyType('char', 'row', {'uniform','gaussian'})
        overrideCoordinateType = symphonyui.core.PropertyType('char', 'row', {'pixels','RF'})
        backgroundIntensity
        overrideRadiiLogical
        meanDisksLogical
        movieMatrix
        xTraj
        yTraj
        timeTraj
        masks
        surroundMask
        specificOpacity
        contrast
        diskExpand
        diskOpacity
        spatialFrequencyPx
        rfSizing
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
            
            % Catch common errors        
            if ~strcmp(obj.trajectory,'natural') && obj.rig.getDevice('Stage').getConfigurationSetting('prerender') == 0
                error('Must have prerender set') % Pre-render required, else trajectory lags and isn't accurate
            end
            
            checkDisks = sum(1:obj.disks);
            checkAssignments = sum([obj.meanDisks obj.contrastDisks obj.meanContrastDisks obj.naturalDisks obj.backgroundDisks]);
            if sum(obj.overrideRadii) > 0
                obj.overrideRadiiLogical = true; % Identifier for later on
                checkDisks = sum(1:size(obj.overrideRadii,2)-1);
                obj.disks = size(obj.overrideRadii,2)-1; % Adjust number of disks.
            end
                       
            if checkDisks > checkAssignments
                error('One or more disks left empty. Please assign a disk type.')
            elseif checkDisks < checkAssignments
                error('Two disk types have been assigned to one location.')
            end
            
            % Empirically defined limits to Stage
            if (obj.disks > 12 && sum([obj.xSliceFrequency obj.ySliceFrequency]) >= 2) ||...
                    (obj.disks > 5 && sum([obj.xSliceFrequency obj.ySliceFrequency]) >= 4) ||...
                    (obj.disks > 3 && sum([obj.xSliceFrequency obj.ySliceFrequency]) >= 6) ||...
                    (obj.disks > 2 && sum([obj.xSliceFrequency obj.ySliceFrequency]) >= 8) 
                if obj.disks == 6 && size(obj.disksIgnoreCut,2) >= 2
                    % Special case that works, can proceed
                else
                    error('Too many disks and slices. Will cause Stage to crash. Please choose a different set.')
                end
            end
            
            % Pull base trajectories and image information.
            [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.imageNo, obj.observerNo,...
                    'amplification', obj.amplification,'mirroring', true);
                
            % Scale pixels in image to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            
            % Produce trajectories
            obj.xTraj = baseMovement.x + fixMovement.x;
            obj.yTraj = baseMovement.y + fixMovement.y;
            
            % We do not need to consider the entire trajectory, however.
            frames = round((obj.stimTime + 50) / 1000 * 200); % max # of frames in DOVES database, with 50ms cushion
            if size(obj.xTraj,2) <= frames
                frames = size(obj.xTraj,2);
            end
            
            % This enables a faster computation time for shorter stimTimes
            obj.xTraj = obj.xTraj(1,1:frames);
            obj.yTraj = obj.yTraj(1,1:frames);
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % DOVES resolution

            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Calculate screen size
            
            if ~strcmp(obj.trajectory,'natural') % Features any disk
                
                % Identify type of disks
                obj.meanDisksLogical = zeros(1,obj.disks);                
                obj.meanDisks(obj.meanDisks == 0) = [];
                obj.meanDisksLogical(obj.meanDisks) = 1;
                      
                % Identify corresponding RF Filter
                [RFFilter,~,obj.rfSizing] = calculateFilter(obj);

                % Prepare trajectory with our RF Filter.
                [wTraj, dist, unwTraj, RFFilterVH] = weightedTrajectory(obj, img2, RFFilter);
                                
                % For an override, convert RF units into pixels.
                if strcmp(obj.overrideCoordinate,'RF') && obj.overrideRadiiLogical == true                   
                    H = obj.rfSizing(2) - obj.rfSizing(1);
                    G = max(canvasSize) / 2 - obj.rfSizing(2);
                    for a = 1:obj.disks+1
                        if obj.overrideRadii(1,a) <= 1
                            obj.overrideRadii(1,a) = obj.overrideRadii(1,a) * obj.rfSizing(1);
                        elseif obj.overrideRadii(1,a) <= 2
                            A = obj.overrideRadii(1,a) - 1; % normalize
                            obj.overrideRadii(1,a) = obj.rfSizing(1) + A * H;
                        elseif obj.overrideRadii(1,a) <= 3
                            A = obj.overrideRadii(1,a) - 2; % normalize
                            obj.overrideRadii(1,a) = obj.rfSizing(2) + A * G;
                        end
                    end
                end

                tic

                % Calculate different disks
                if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                    movieOutput = linearEquivalencyProjectionMovie(obj, wTraj, dist, RFFilterVH, unwTraj);
                else
                    movieOutput = linearEquivalencySlicedProjectionMovie(obj, wTraj, dist, RFFilterVH, unwTraj);
                end
                
                toc
                
                obj.movieMatrix = uint8(cat(4,unwTraj,movieOutput));
            end
            
            % There may be leaky pixels around the edge that could make
            % comparisons difficult. So, we build an mask to keep comparison controlled.
            [xx, yy] = meshgrid(1:2*canvasSize(1),1:2*canvasSize(2));
            m = sqrt((xx-canvasSize(1)).^2+(yy-canvasSize(2)).^2);
            obj.surroundMask = m >= (max(obj.radii).*(1-obj.diskExpand/100));
            protectiveMask = zeros(2*canvasSize(2),2*canvasSize(1));
            protectiveMask(round(canvasSize(2) - ceil(canvasSize(2)/2) + 1) : round(canvasSize(2) + ceil(canvasSize(2)/2) - 1), ...
                round(canvasSize(1) - ceil(canvasSize(1)/2)) + 1 : round(canvasSize(1) + ceil(canvasSize(1)/2)) - 1) = 1;
            protectiveMask = abs(protectiveMask - 1);
            obj.surroundMask = obj.surroundMask + protectiveMask; 
            
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
            epoch.addParameter('radii', obj.radii); % in pixels
            epoch.addParameter('rfSize',obj.rfSizing); % in pixels
            
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
            scene = stage.builtin.stimuli.Movie(obj.movieMatrix);
            scene.size = [size(obj.movieMatrix,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.movieMatrix,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            
            % Apply eye trajectories to move image around.
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            p.addStimulus(scene);

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)rem(state.time,cycleTime) >= obj.preTime * 1e-3 && rem(state.time,cycleTime) < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            % Apply far-surround mask
            surround = stage.builtin.stimuli.Rectangle();
            surround.position = canvasSize/2;
            surround.size = [2*canvasSize(1) 2*canvasSize(2)];
            surround.color = obj.backgroundIntensity;
            annulusS = uint8(obj.surroundMask*255);
            surroundMaskX = stage.core.Mask(annulusS);
            surround.setMask(surroundMaskX);
            p.addStimulus(surround);
            
        end
        
        % Calculate RF Filter
        function [g, b, f] = calculateFilter(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size

            % Convert to pixels
            centerSigmaPix = obj.rfSigmaCenter / obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            surroundSigmaPix = obj.rfSigmaSurround / obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');

            centerGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],centerSigmaPix);
            surroundGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],surroundSigmaPix);

            % Calculate difference of gaussians
            diffGaussian.raw = centerGaus - surroundGaus;

            % We want to make our background light level > 0 s.t. we can
            % observe inhibitory effects
            filte = diffGaussian.raw ./ max(diffGaussian.raw(:)); % pre-normalize filter
            supp = abs(min(filte(:))); % pre-normalize suppressive region

            filte = filte+supp; % scaled  s.t. no zeros exist

            g = filte./max(filte(:)); % re-normalize filter
            b = supp / max(filte(:)); % re-normalize background intensity
            
            % We can use this moment to calculate the size of our RFs.
            threshold = b ./ 5; % threshold by which we define our surround.

            slice = g(round(median(1:size(g,1))),round(median(1:size(g,2))):end); % take a 2D slice

            sizing = slice <= threshold; % find annulus

            len = 1:size(sizing,2);
            centerRes = sizing .* len; % make values index
            centerRes(centerRes == 0) = [];
            
            % This parameter determines the sizing of the center/surround RF
            f(1,1) = min(centerRes); % in pixels
            f(1,2) = max(centerRes);
        end
        
        % Apply RF Filter over the entire trajectory.
        function [r, m, u, o] = weightedTrajectory(obj, img, RFFilter)
            
            % calculate trajectory size
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            imgSize = ceil(canvasSize / (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')));
            xRange = floor(imgSize(1) / 2);
            yRange = floor(imgSize(2) / 2);
            
            % resize filter accordingly
            tempFilt = RFFilter(round(size(RFFilter,1)/2) - yRange : round(size(RFFilter,1)/2) + yRange,...
                round(size(RFFilter,2)/2) - xRange : round(size(RFFilter,2)/2) + xRange);
            
            % create mesh grid
            [xx, yy] = meshgrid(1:imgSize(1),1:imgSize(2));
            m = sqrt((xx-imgSize(1)/2).^2+(yy-imgSize(2)/2).^2);
            
            % check that m and tempFilt are same size
            if size(m) ~= size(tempFilt)
                tempFilt = imresize(tempFilt,size(m));
            end
            
            % while the rig automatically centers the stimulus, our
            % calculation doesn't.
            centering = obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'); % in mu
            centeringPix = centering ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            centeredXTraj = round(obj.xTraj - centeringPix(1));
            centeredYTraj = round(obj.yTraj + centeringPix(2));
            
            % create trajectory  
            r = zeros(size(m,1),size(m,2),1,size(obj.xTraj,2));
            u = zeros(size(m,1),size(m,2),1,size(obj.xTraj,2));
            for a = 1:size(obj.xTraj,2)
                % create image
               u(:,:,1,a) = img(centeredYTraj(1,a)-yRange:centeredYTraj(1,a)+yRange-1,...
                   centeredXTraj(1,a)-xRange:centeredXTraj(1,a)+xRange-1); 
               r(:,:,1,a) = u(:,:,1,a) .* tempFilt;
            end
            
            o = tempFilt;
        end
        
        % Calculate linear equivalency without slices
        function outputMovie = linearEquivalencyProjectionMovie(obj, r, m, RFFilterVH, unwTraj)
            
            outputMovie = unwTraj; % typical trajectory
            imgSize = [size(r,1) size(r,2)];
            minimumPixels = 12; % minimum amount of pixels for us to get a reasonable statistic.

            radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            
            if obj.overrideRadiiLogical == true
                radius = obj.overrideRadii ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % in VH pixels
            end
                       
            if sum(obj.meanDisksLogical) > 0 % Check whether the whole trajectory is needed.
                % Calculate statistics
                for a = 1:size(radius,2) - 1
                    if obj.meanDisksLogical(1,a) > 0 % Only regions where we absolutely need to calculate (speed increase)
                        filt = m >= radius(1,a) & m <= radius(1,a+1); % Logical mask.
                        
                        % check enough pixels fill the space
                        check = filt;
                        check(check == 0) = [];
                        if size(check,2) < minimumPixels % Not enough pixels to continue
                            error('Not enough pixels in region. Please change mask radii.')
                        end
                        
                        regionSize = sum(filt(:)); % pixels in area
                        
                        if strcmp(obj.meanIntegration,'gaussian')
                            T = RFFilterVH .* filt;
                            R = sum(T(:)) / regionSize * 255; % Mean equivalence
                        end

                        % Apply filter across full trajectory.
                        for n = 1:size(r,4)
                            
                            if obj.meanDisksLogical(1,a) == 1 % Disk belongs
                                
                                if strcmp(obj.meanIntegration,'gaussian')
                                    S = r(:,:,1,n) .* filt; % RF weighted trajectory
                                else
                                    S = unwTraj(:,:,1,n) .* filt; % Unweighted trajectory
                                end
                                
                                if isempty(S) % Becomes an issue when pixel weight is very small (approx 0).
                                    S = 0;
                                end
                                
                                if strcmp(obj.meanIntegration,'gaussian')
                                    equivalency = sum(S(:)) / (regionSize  * R); % Renormalize
                                else
                                    equivalency = sum(S(:)) / (regionSize * 255); % Raw average
                                end
                                
                                deltaDisk = (equivalency - (obj.backgroundIntensity * 255)) .* filt; % change relative to global mean
                                outputMovie(:,:,1,n) = outputMovie(:,:,1,n) + deltaDisk; % add filter
                            end
                        end
                    end
                end
            end
        end
        
        % Calculate linear equivalency with slices
        function outputMovie = linearEquivalencySlicedProjectionMovie(obj, r, m, RFFilterVH, unwTraj)
            
            outputMovie = unwTraj; % typical trajectory
            imgSize = [size(r,1) size(r,2)];
            minimumPixels = 12; % minimum amount of pixels for us to get a reasonable statistic.
            
            [xx,yy] = meshgrid(1:imgSize(2),1:imgSize(1));
            k = atan((xx - imgSize(2)/2) ./ (yy - imgSize(1)/2)); % Calculate theta
            k = k ./ max(k(:)) .* 90; % Convert to degrees
            k = abs(k - 90); % Rotate for proper polar coordinates
            k(1:floor(imgSize(1)/2),:) = k(1:floor(imgSize(1)/2),:) + 180;
            
            k = mod(k - obj.rotateSlices,360); % rotate as needed.

            radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            
            if obj.overrideRadiiLogical == true
                radius = obj.overrideRadii ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % in VH pixels
            end
            
            % Calculate angles to operate over
            if obj.xSliceFrequency > 0 && obj.ySliceFrequency == 0
                thetaX = 90 / obj.xSliceFrequency;
                theta = 0:thetaX:90-(1E-10);
                theta = [theta theta + 180 360];
            elseif obj.xSliceFrequency == 0 && obj.ySliceFrequency > 0
                thetaY = 90 / obj.ySliceFrequency;
                theta = 90:thetaY:180-(1E-10);
                theta = [theta theta + 180 450];
                k(round(imgSize(1)/2):end,round(imgSize(2)/2):end) = k(round(imgSize(1)/2):end,round(imgSize(2)/2):end) + 360; % add
            elseif obj.xSliceFrequency > 0 && obj.ySliceFrequency > 0
                thetaX = 90 / obj.xSliceFrequency;
                thetaY = 90 / obj.ySliceFrequency;
                theta = [0:thetaX:90-(1E-10) 90:thetaY:180-(1E-10)];
                theta = [theta theta + 180 360];
            end
                        
            if sum(obj.meanDisksLogical) > 0 % Check whether the whole trajectory is needed.
                % Calculate statistics
                for a = 1:size(radius,2) - 1
                    for b = 1:size(theta,2) - 1
                        if obj.meanDisksLogical(1,a) > 0 % Only regions where we absolutely need to calculate (speed increase)
                            radiusFilt = m >= radius(1,a) & m <= radius(1,a+1); % Radial filter (r)
                            angFilt = k >= theta(1,b) & k <= theta(1,b+1); % Angular filter (theta)

                            if ismember(a,obj.disksIgnoreCut) % Ignore angular filter
                                angFilt = ones(size(angFilt));
                            end  

                            filt = radiusFilt .* angFilt;
                            check = filt .* m;
                            check(check == 0) = [];

                            if size(check,2) < minimumPixels % Not enough pixels to continue
                                error('Not enough pixels in region. Please change mask radii.')
                            end
                            
                            regionSize = sum(filt(:));
                            
                            if strcmp(obj.meanIntegration,'gaussian')
                                T = RFFilterVH .* filt;
                                R = sum(T(:)) / regionSize * 255; % Mean equivalence
                            end

                            % Apply filter across full trajectory.
                            for n = 1:size(r,4)

                                if obj.meanDisksLogical(1,a) == 1 % Disk belongs

                                    if strcmp(obj.meanIntegration,'gaussian')
                                        S = r(:,:,1,n) .* filt; % RF weighted trajectory
                                    else
                                        S = unwTraj(:,:,1,n) .* filt; % Unweighted trajectory
                                    end

                                    if isempty(S) % Becomes an issue when pixel weight is very small (approx 0).
                                        S = 0;
                                    end

                                    if strcmp(obj.meanIntegration,'gaussian')
                                        equivalency = sum(S(:)) / (regionSize  * R); % Renormalize
                                    else
                                        equivalency = sum(S(:)) / (regionSize  * 255); % Raw average
                                    end
                                    
                                    deltaDisk = (equivalency - (obj.backgroundIntensity * 255)) .* filt; % change relative to global mean
                                    outputMovie(:,:,1,n) = outputMovie(:,:,1,n) + deltaDisk; % add filter

                                end

                            end
                        end
                    end
                end
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