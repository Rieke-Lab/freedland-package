% Replace a natural movie with a variety of linear equivalent disks.
classdef RFDiskArray < edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol
    properties
        % stimulus times
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % natural image trajectory
        imageNo = 1; % natural image number (1 to 11)
        observerNo = 1; % observer number (1 to 19)
        amplification = 1; % amplify fixations by X.  
        trajectory = 'natural'; % which type of stimulus to present: natural image, linear disk replacement ,or both?        
        
        % RF field information
        rfSigmaCenter = 60; % (um) enter from difference of gaussians fit, for overlaying RF
        rfSigmaSurround = 160; % (um) enter from difference of gaussians fit, for overlaying RF
        
        % disk information
        disks = 5; % number of disks that replace the natural image.
        diskFocus = 'none'; % where to place the most masks.
        xSliceFrequency = 1; % how many slices to cut between 0 and 90 degrees.
        ySliceFrequency = 1; % how many slices to cut between 90 and 180 degrees.
        disksIgnoreCut = [1 2]; % starting from the center disk and moving outwards, how many disks should we 'forget' to slice?
        disksIgnoreMask = [0 0]; % starting from the center disk and moving outwards, how many disks should we entirely 'forget' to place? (reveals natural image underneath).
        diskEvenness = 1; % alters width of disks. choose 0.1 - 10. 10 = even widths, 0.1 = very uneven widths. Must be >0.
        minimumPixels = 12; % minimum number of pixels required to calculate average. if under, we will extend our disk width to reach this value.
        diskExpand = 0.5; % percent to increase each disk's radius (slight overlaps ensure masks cover natural image completely)
        
        % additional parameters
        offsetWidth = 0 % in microns. offsets trajectory (x direction) to expose new areas of focus.
        offsetHeight = 0 % in microns. offsets trajectory (y direction) to expose new areas of focus.
        mirroring = true; % mirror image when reaching an edge.
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        diskFocusType = symphonyui.core.PropertyType('char', 'row', {'none','center','near-center','center/near-surround','near-surround','near/far surround','far-surround'})
        trajectoryType = symphonyui.core.PropertyType('char', 'row', {'natural','disk','both'})
        backgroundIntensity
        imageMatrix
        xTraj
        yTraj
        timeTraj
        linear
        radii
        masks
        surroundMask
        theta
        specificOpacity
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)
                   
            % For Symphony
            prepareRun@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj);
            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            if strcmp(obj.trajectory,'both')
                obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',true);
            else
                obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            end
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % ensure re-render is set.
            if ~strcmp(obj.trajectory,'natural') && obj.rig.getDevice('Stage').getConfigurationSetting('prerender') == 0
                error('Must have prerender set')
            end
            
            % Identify image
            imageIdentifier = [5 8 13 17 23 27 31 39 56 64 100];
            imageVal = imageIdentifier(obj.imageNo);
            
            % Convert um to pixels
            obj.offsetHeight = obj.rig.getDevice('Stage').um2pix(obj.offsetHeight);
            obj.offsetWidth = obj.rig.getDevice('Stage').um2pix(obj.offsetWidth);
            
            % Make paths for every amplification.
            [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal, obj.observerNo,...
                    'amplification', obj.amplification,'offSetHeight', obj.offsetHeight,'offSetWidth',...
                    obj.offsetWidth, 'mirroring', obj.mirroring);
                
            % Adjust image
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix = uint8(img2);
            
            % Recombine trajectories
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

            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % calculate screen size
            
            if ~strcmp(obj.trajectory,'natural') % features any linear equivalency.
                % Identify corresponding RF Filter
                [RFFilter,~,sizing] = calculateFilter(obj);

                % Prepare trajectory with our RF Filter.
                [traj, dist] = weightedTrajectory(obj, img2, RFFilter);

                % Calculate linear equivalency
                if obj.xSliceFrequency == 0 && obj.ySliceFrequency == 0
                    [obj.linear, obj.radii] = linearEquivalency(obj, traj, dist, sizing);
                else
                    [obj.linear, obj.radii, obj.theta] = linearEquivalencySliced(obj, traj, dist, sizing);
                end

                obj.radii(isnan(obj.radii)) = []; % remove insufficient radii

                % Make masks
                obj.radii = obj.radii .* (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % convert from VH to px
                [xx, yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
                m = sqrt((xx-canvasSize(1)/2).^2+(yy-canvasSize(2)/2).^2); % recalculate for mask size
                k = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2)); % calculate theta
                k = k ./ max(k(:)) .* 90; % convert to degrees
                k = abs(k - 90); % rotate for proper coordinates
                k(1:round(canvasSize(2)/2)-1,:) = k(1:round(canvasSize(2)/2)-1,:) + 180; % finish off
                
                if obj.xSliceFrequency == 0 && obj.ySliceFrequency > 0 % special case
                    k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) = k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) + 360; % add
                end

                % The case with no slices
                if obj.xSliceFrequency == 0 && obj.ySliceFrequency == 0
                    obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2));
                    obj.specificOpacity = ones(1,size(obj.masks,3));
                    for a = 1:size(obj.radii,2) - 1
                        obj.masks(:,:,a) = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                            & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                    end
                    
                    for so = 1:size(obj.disksIgnoreMask,2)
                        r = obj.disksIgnoreMask(1,so);
                        if r > 0 % check real index.
                            obj.specificOpacity(1,r) = 0; % set ignored masks to opacity 0.
                        end
                    end
                    
                else % The case with slices
                    obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2),size(obj.theta,2) - 1);
                    obj.specificOpacity = ones(1,size(obj.masks,3),size(obj.masks,4));
                    for a = 1:size(obj.radii,2) - 1
                        for b = 1:size(obj.theta,2) - 1
                            dist = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                                & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                            th = k >= (obj.theta(1,b)) & k <= (obj.theta(1,b+1));
       
                            % ignore cuts in the center
                            if ismember(a,obj.disksIgnoreCut);
                                th = ones(size(k));
                            end
                            
                            obj.masks(:,:,a,b) = dist .* th;
                        end
                    end
                    
                    for so = 1:size(obj.disksIgnoreMask,2)
                        r = obj.disksIgnoreMask(1,so);
                        if r > 0 % check real index.
                            obj.specificOpacity(1,r,:) = 0; % set ignored masks to opacity 0.
                        end
                    end
                end                
            else
                obj.radii = max(canvasSize) / 2;
            end
            
            % Build surround mask to ensure no natural image pixels leak around edges.
            [xx, yy] = meshgrid(1:2*canvasSize(1),1:2*canvasSize(2));
            m = sqrt((xx-canvasSize(1)).^2+(yy-canvasSize(2)).^2);
            obj.surroundMask = m >= (max(obj.radii).*(1-obj.diskExpand/100));
            protectiveMask = zeros(2*canvasSize(2),2*canvasSize(1));
            protectiveMask(round(canvasSize(2)) - ceil(canvasSize(2)/2) : round(canvasSize(2)) + ceil(canvasSize(2)/2), ...
                round(canvasSize(1)) - ceil(canvasSize(1)/2) : round(canvasSize(1)) + ceil(canvasSize(1)/2)) = 1;
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
                duration = 2* (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            end
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('radii', obj.radii);
            
            % we will also add some metadata from Stage.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
        end
        
        function p = createPresentation(obj)
            % Prep stage for presentation
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
                                    
            if ~strcmp(obj.trajectory,'both')
                p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            else % present both stimuli in succession
                p = stage.core.Presentation(2 * (obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            end

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
            
            % Apply eye trajectories to move image around.
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, rem(state.time - obj.preTime/1e3,cycleTime), p0));
            
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
            p.addController(scenePosition);

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)rem(state.time,cycleTime) >= obj.preTime * 1e-3 && rem(state.time,cycleTime) < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);

            %%%%%% Apply masks %%%%%%
            
            % always apply far-surround mask (avoid leaking pixels on edge)
            surround = stage.builtin.stimuli.Rectangle();
            surround.position = canvasSize/2;
            surround.size = [2*canvasSize(1) 2*canvasSize(2)];
            surround.color = obj.backgroundIntensity;
            annulusS = uint8(obj.surroundMask*255);
            surroundMaskX = stage.core.Mask(annulusS);
            surround.setMask(surroundMaskX);
            p.addStimulus(surround);
            
            if ~strcmp(obj.trajectory,'natural') 
                
                if obj.xSliceFrequency == 0 && obj.ySliceFrequency == 0

                    % no slices, assumes radial symmetry.
                    for q = 1:size(obj.linear,1)
                        annulus = stage.builtin.stimuli.Rectangle();
                        annulus.position = canvasSize/2;
                        annulus.size = [canvasSize(1) canvasSize(2)];
                        annulus.opacity = obj.specificOpacity(1,q);
                        annulusA = uint8(obj.masks(:,:,q)*255);
                        annulusMaskX = stage.core.Mask(annulusA);
                        annulus.setMask(annulusMaskX);

                        if ~strcmp(obj.trajectory,'both') 
                            annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, obj.linear(q,:)));
                        else
                            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
                            annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, obj.linear(q,:)));
                            
                            sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                    @(state)state.time >= cycleTime);
                            p.addController(sceneVisible);
                        end

                        p.addStimulus(annulus);
                        p.addController(annulusLED);              
                    end
                else
                    % with slices, more complex averaging.
                    for q = 1:size(obj.linear,1)
                        for s = 1:size(obj.linear,2)
                            annulus = stage.builtin.stimuli.Rectangle();
                            annulus.position = canvasSize/2;
                            annulus.size = [canvasSize(1) canvasSize(2)];
                            annulus.opacity = obj.specificOpacity(1,q,s);
                            annulusA = uint8(obj.masks(:,:,q,s)*255);
                            annulusMaskX = stage.core.Mask(annulusA);
                            annulus.setMask(annulusMaskX);
                            
                            F = obj.linear(q,s,:);
                            F = reshape(F, [1 size(F,3)]);

                            if ~strcmp(obj.trajectory,'both') 
                                annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, F));
                            else
                                cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
                                annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, F));
                                
                                sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                    @(state)state.time >= cycleTime);
                                p.addController(sceneVisible);
                            end

                            p.addStimulus(annulus);
                            p.addController(annulusLED);     
                        end
                    end
                end
            end
            
            % Match mask color to linear equivalency in real time.
            function s = getBackground(obj, time, equivalency)
                if time < 0
                    s = obj.backgroundIntensity;
                elseif time > obj.timeTraj(end)
                    s = equivalency(1,end); % hold last color
                else
                    s = interp1(obj.timeTraj,equivalency,time);
                end
            end
        end
        
        % Calculate RF Filter
        function [g, b, f] = calculateFilter(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size

            % Convert to pixels
            centerSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaCenter);
            surroundSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaSurround);

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
            f(1,1) = min(centerRes);
            f(1,2) = max(centerRes);
        end
        
        % Apply RF Filter over the entire trajectory.
        function [r, m] = weightedTrajectory(obj, img, RFFilter)
            
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
            
            % create trajectory  
            r = zeros(size(m,1),size(m,2),1,size(obj.xTraj,2));
            for a = 1:size(obj.xTraj,2)
                % create image
               tempImg = img(obj.yTraj(1,a)-yRange:obj.yTraj(1,a)+yRange,obj.xTraj(1,a)-xRange:obj.xTraj(1,a)+xRange); 
               r(:,:,1,a) = tempImg .* tempFilt;
            end
        end
        
        % Calculate linear equivalency without slices
        function [q, radius] = linearEquivalency(obj, r, m, sizing)
            
            q = zeros(obj.disks,size(r,4));
            imgSize = [size(r,1) size(r,2)];
            
            % identify center, near-, and far-surround
            if strcmp(obj.diskFocus,'center')
                cent =  sizing(1,1) / 3;
            elseif strcmp(obj.diskFocus,'near-center')
                cent = sizing(1,1) / 2;
            elseif strcmp(obj.diskFocus,'center/near-surround')
                cent = sizing(1,1);
            elseif strcmp(obj.diskFocus,'near-surround')
                cent = ((sizing(1,2) - sizing(1,1)) / 2) + sizing(1,1);
            elseif strcmp(obj.diskFocus,'near/far surround')
                cent = sizing(1,2);
            elseif strcmp(obj.diskFocus,'far-surround')
                cent = ((max(m) - sizing(1,2)) / 2) + sizing(1,2);
            end
            
            % calculate distribution of masks to calculate mean over
            if ~strcmp(obj.diskFocus,'none')
                maskDistribution = normpdf(1:max(imgSize)/(obj.disks):max(imgSize),cent,cent ./ obj.diskEvenness);
                maskDistribution = maskDistribution ./ (max(maskDistribution(:))) - (min(maskDistribution(:))); % normalize
                maskDistribution = abs(maskDistribution - 1);
                maskDistribution = maskDistribution ./ sum(maskDistribution(:)) .* max(imgSize);
                radius = [0 cumsum(maskDistribution)] ./ 2; % this gives us the radius of each mask
            else
                radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            end
            
            % calculate mean
            for a = 1:size(radius,2) - 1
                filt = m >= radius(1,a) & m <= radius(1,a+1); % calculate filter based on radius
                check = filt .* m;
                check(check == 0) = [];
                
                if size(check,2) <= obj.minimumPixels % must have enough pixels!
                    temp = m >= radius(1,a); % expand outer edge
                    temp = temp .* m;
                    diskDiff = temp - (filt .* m); % subtract the larger disk from the smaller
                    diskDiff(diskDiff == 0) = [];
                    diskDiff2 = sort(diskDiff);
                    
                    if size(diskDiff2,2) >= obj.minimumPixels
                        radius(1,a+1) = diskDiff2(1,obj.minimumPixels); % find the smallest radius that abides by our pixel minimum
                        filt = m >= radius(1,a) & m <= radius(1,a+1); % recalculate filter
                    else % at edge of image
                        radius(1,a+1) = max(imgSize);
                        filt = zeros(size(r,1),size(r,2)); % will cause the mean to be the same as background
                    end
                        
                end
                
                % use filter across entire trajectory.
                for n = 1:size(r,4)
                    S = r(:,:,1,n) .* filt;
                    S(S==0) = [];
                    
                    if isempty(S)
                        S = 0;
                    end
                    
                    q(a,n) = mean(S) / 255;
                end
            end
        end
        
        % Calculate linear equivalency with slices
        function [q, radius, theta] = linearEquivalencySliced(obj, r, m, sizing)
            
            imgSize = [size(r,1) size(r,2)];
            
            [xx,yy] = meshgrid(1:imgSize(2),1:imgSize(1));
            k = atan((xx - imgSize(2)/2) ./ (yy - imgSize(1)/2)); % calculate theta
            k = k ./ max(k(:)) .* 90; % convert to degrees
            k = abs(k - 90); % rotate for proper coordinates
            k(1:floor(imgSize(1)/2),:) = k(1:floor(imgSize(1)/2),:) + 180; % finish off
            
            % identify center, near-, and far-surround
            if strcmp(obj.diskFocus,'center')
                cent =  sizing(1,1) / 3;
            elseif strcmp(obj.diskFocus,'near-center')
                cent = sizing(1,1) / 2;
            elseif strcmp(obj.diskFocus,'center/near-surround')
                cent = sizing(1,1);
            elseif strcmp(obj.diskFocus,'near-surround')
                cent = ((sizing(1,2) - sizing(1,1)) / 2) + sizing(1,1);
            elseif strcmp(obj.diskFocus,'near/far surround')
                cent = sizing(1,2);
            elseif strcmp(obj.diskFocus,'far-surround')
                cent = ((max(m) - sizing(1,2)) / 2) + sizing(1,2);
            end
            
            % calculate distribution of masks to calculate mean over
            if ~strcmp(obj.diskFocus,'none')
                maskDistribution = normpdf(1:max(imgSize)/(obj.disks):max(imgSize),cent,cent ./ obj.diskEvenness);
                maskDistribution = maskDistribution ./ (max(maskDistribution(:))) - (min(maskDistribution(:))); % normalize
                maskDistribution = abs(maskDistribution - 1);
                maskDistribution = maskDistribution ./ sum(maskDistribution(:)) .* max(imgSize);
                radius = [0 cumsum(maskDistribution)] ./ 2; % this gives us the radius of each mask
            else
                radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            end
            
            % calculate angles to operate over
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
            
            q = zeros(obj.disks,size(theta,2)-1,size(r,4)); % in form [disk #, angle, mean value]
            
            % calculate mean
            for a = 1:size(radius,2) - 1
                for b = 1:size(theta,2) - 1
                    radiusFilt = m >= radius(1,a) & m <= radius(1,a+1); % calculate filter based on radius
                    angFilt = k >= theta(1,b) & k <= theta(1,b+1);
                    
                    if ismember(a,obj.disksIgnoreCut)
                        angFilt = ones(size(angFilt));
                    end  
                    
                    filt = radiusFilt .* angFilt;
                    check = filt .* m;
                    check(check == 0) = [];

                    if size(check,2) <= obj.minimumPixels % must have enough pixels!
                        temp = m >= radius(1,a); % expand outer edge
                        temp = temp .* angFilt .* m;
                        diskDiff = temp - (filt .* m); % subtract the larger disk from the smaller
                        diskDiff(diskDiff == 0) = [];
                        diskDiff2 = sort(diskDiff);

                        if size(diskDiff2,2) >= obj.minimumPixels % check we're not at the end of the image
                            radius(1,a+1) = diskDiff2(1,obj.minimumPixels); % find the smallest radius that abides by our pixel minimum
                            radiusFilt = m >= radius(1,a) & m <= radius(1,a+1); % recalculate filter
                            filt = radiusFilt .* angFilt;
                        else % at edge of image
                            radius(1,a+1) = NaN; % insufficient radius
                            filt = zeros(size(r,1),size(r,2)); % mean = 0
                        end

                    end

                    % use filter across entire trajectory.
                    for n = 1:size(r,4)
                        S = r(:,:,1,n) .* filt;
                        S(S==0) = [];
                        
                        if isempty(S)
                            S = 0;
                        end
                        
                        q(a,b,n) = mean(S) / 255;
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