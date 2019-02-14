% Amplify fixations with a carefully probed visual field.
classdef RFDiskArray < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        imageNo = 1; % natural image number (1 to 11)
        observerNo = 1; % observer number (1 to 19)
        amplification = 1; % amplify fixations by X.        
        rfSigmaCenter = 60; % (um) enter from difference of gaussians fit, for overlaying RF
        rfSigmaSurround = 160; % (um) enter from difference of gaussians fit, for overlaying RF 
        disks = 10; % number of disks in RF.
        diskFocus = 'none'; % where to place the most masks.
        xAxisSlice = false; % slice all our masks horizontally. Consider each section independently.
        yAxisSlice = false; % slice all our masks vertically. Consider each section independently.
        diskEvenness = 1; % Choose 0.1 - 10. 10 = even widths, 0.1 = very uneven widths. Must be >0.
        minimumPixels = 8; % minimum number of pixels required to calculate average.
        diskExpand = 0.5; % percent to increase each disk's radius (slight overlaps ensure masks cover natural image completely)
        boost = false; % apply a light level boost to a region?
        boostRegion = [100 300]; % boost levels between both radii (in um).
        boostRegionBy = 0.2; % add (between 0 and 1) light intensity to region.
         
        offsetWidth = 0 % in microns
        offsetHeight = 0 % in microns
        
        mirroring = true; % mirror image when reaching an edge.
        
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(1) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        diskFocusType = symphonyui.core.PropertyType('char', 'row', {'none','center','near-center','center/near-surround','near-surround','near/far surround','far-surround'})
        backgroundIntensity
        imageMatrix
        xTraj
        yTraj
        timeTraj
        linear
        radii
        masks
        surroundMask
        counter
        theta
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
            % to help refine parameters.
            obj.xTraj = obj.xTraj(1,1:frames);
            obj.yTraj = obj.yTraj(1,1:frames);
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % DOVES resolution
              
            % Identify corresponding RF Filter
            [RFFilter,~,sizing] = calculateFilter(obj);
            
            % Prepare trajectory with our RF Filter.
            [traj, dist] = weightedTrajectory(obj, img2, RFFilter);
            
            % Calculate linear equivalency
            if obj.xAxisSlice == false && obj.yAxisSlice == false
                [obj.linear, obj.radii] = linearEquivalency(obj, traj, dist, sizing);
            else
                [obj.linear, obj.radii, obj.theta] = linearEquivalencySliced(obj, traj, dist, sizing);
            end
            
            obj.radii(isnan(obj.radii)) = []; % remove insufficient radii
            
            % Make masks
            obj.radii = obj.radii .* (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % convert from VH to px
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); 
            [xx, yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            m = sqrt((xx-canvasSize(1)/2).^2+(yy-canvasSize(2)/2).^2); % recalculate for mask size
            k = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2)); % calculate theta
            k = k ./ max(k(:)) .* 90; % convert to degrees
            k = abs(k - 90); % rotate for proper coordinates
            k(1:round(canvasSize(2)/2)-1,:) = k(1:round(canvasSize(2)/2)-1,:) + 180; % finish off
            
            if obj.xAxisSlice == false && obj.yAxisSlice == true
                k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) = k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) + 360; % add
            end
           
            % The case with no slices
            if obj.xAxisSlice == false && obj.yAxisSlice == false
                obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2));
                for a = 1:size(obj.radii,2) - 1
                    obj.masks(:,:,a) = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                        & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                    % Apply boost if desired
                    if obj.boost == true
                        if obj.radii(1,a) >= obj.rig.getDevice('Stage').um2pix(obj.boostRegion(1,1))&& ...
                                obj.radii(1,a+1) <= obj.rig.getDevice('Stage').um2pix(obj.boostRegion(1,2))
                            obj.linear(a,:) = obj.linear(a,:) + obj.boostRegionBy;
                        end
                    end
                end
            else % The case with slices
                obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2),size(obj.theta,2) - 1);
                for a = 1:size(obj.radii,2) - 1
                    for b = 1:size(obj.theta,2) - 1
                        dist = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                            & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                        th = k >= (obj.theta(1,b)) & k <= (obj.theta(1,b+1));
                        obj.masks(:,:,a,b) = dist .* th;
                        % apply boost if desired
                        if obj.boost == true
                            if obj.radii(1,a) >= obj.rig.getDevice('Stage').um2pix(obj.boostRegion(1,1))&& ...
                                    obj.radii(1,a+1) <= obj.rig.getDevice('Stage').um2pix(obj.boostRegion(1,2))
                                obj.linear(a,b,:) = obj.linear(a,b,:) + obj.boostRegionBy;
                            end
                        end
                    end
                end
            end
            
            % Build surround mask to ensure no natural image pixels leak around edges.
            [xx, yy] = meshgrid(1:2*canvasSize(1),1:2*canvasSize(2));
            m = sqrt((xx-canvasSize(1)).^2+(yy-canvasSize(2)).^2);
            obj.surroundMask = m >= (max(obj.radii).*(1-obj.diskExpand/100));
            
            % Adjust axes and units for monitor (VH Pixels)
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            obj.xTraj = obj.xTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj = obj.yTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
        
            % Only apply masks on second run-thru
            obj.counter = 0;
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
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
            
            % Apply eye trajectories to move image around
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3, p0));
            
            function p = getScenePosition(obj, time, p0)
                if time < 0
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
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);

            %%%%%% Apply masks %%%%%%
            
            if obj.counter == 1 % second-run through
            
                % background mask (avoid leaking pixels on edge)
                surround = stage.builtin.stimuli.Rectangle();
                surround.position = canvasSize/2;
                surround.size = [2*canvasSize(1) 2*canvasSize(2)];
                surround.color = obj.backgroundIntensity;
                annulusS = uint8(obj.surroundMask*255);
                surroundMaskX = stage.core.Mask(annulusS);
                surround.setMask(surroundMaskX);
                p.addStimulus(surround);
                
                if obj.xAxisSlice == false && obj.yAxisSlice == false

                    % no slices, assumes radial symmetry.
                    for q = 1:size(obj.linear,1)
                        annulus = stage.builtin.stimuli.Rectangle();
                        annulus.position = canvasSize/2;
                        annulus.size = [canvasSize(1) canvasSize(2)];
                        annulusA = uint8(obj.masks(:,:,q)*255);
                        annulusMaskX = stage.core.Mask(annulusA);
                        annulus.setMask(annulusMaskX);

                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                            'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, obj.linear(q,:)));

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
                            annulusA = uint8(obj.masks(:,:,q,s)*255);
                            annulusMaskX = stage.core.Mask(annulusA);
                            annulus.setMask(annulusMaskX);
                            
                            F = obj.linear(q,s,:);
                            F = reshape(F, [1 size(F,3)]);

                            annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, F));

                            p.addStimulus(annulus);
                            p.addController(annulusLED);     
                        end
                    end
                end
            end
            
            % Match mask color to linear equivalency in real time.
            function s = getBackground(obj, time, equivalency)
                if time < 0 || time > obj.timeTraj(end)
                    s = obj.backgroundIntensity;
                else 
                    s = interp1(obj.timeTraj,equivalency,time);
                end
            end
            
            obj.counter = mod(obj.counter + 1,2);
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
            
            if obj.xAxisSlice == true && obj.yAxisSlice == false
                theta = [0 180 360];
            elseif obj.xAxisSlice == false && obj.yAxisSlice == true
                theta = [90 270 450];
                k(round(imgSize(2)/2):end,round(imgSize(1)/2):end) = k(round(imgSize(2)/2):end,round(imgSize(1)/2):end) + 360; % add
            else
                theta = [0 90 180 270 360];
            end
            
            q = zeros(obj.disks,size(theta,2)-1,size(r,4)); % in form [disk #, angle, mean value]
            
            % calculate mean
            for a = 1:size(radius,2) - 1
                for b = 1:size(theta,2) - 1
                    radiusFilt = m >= radius(1,a) & m <= radius(1,a+1); % calculate filter based on radius
                    angFilt = k >= theta(1,b) & k <= theta(1,b+1);
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
            tf = obj.numEpochsPrepared <= obj.numberOfAverages * 2 - 1;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted <= obj.numberOfAverages * 2 - 1;
        end
    end
end