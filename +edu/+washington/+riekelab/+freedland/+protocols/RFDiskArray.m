% Amplify fixations with a carefully probed visual field.
classdef RFDiskArray < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        imageNo = 1; % natural image number (1 to 11)
        observerNo = 1; % observer number (1 to 19)
        amplification = 1; % amplify fixations by X.
        
        rfSigmaCenter = 60; % (um) enter from difference of gaussians fit, uses to calculate center size
        rfSigmaSurround = 160; % (um) enter from difference of gaussians fit, uses to calculate center size
        disks = 20; % number of disks in RF. Must be a whole number.
        diskFocus = 'none'; % where to place the most masks
        diskEvenness = 1; % Choose 0.1 - 10. 10 = even widths, 0.1 = very uneven widths. Must be >0
        diskExpand = 0.5; % percent to increase each disk's radius (slight overlaps ensure masks cover natural image completely)
        boost = false; % apply a light level boost to a region?
        boostRegion = [100 300]; % boost levels between both radii (in um).
        boostRegionBy = 0.2; % add (between 0 and 1) light intensity to region.
         
        offsetWidth = 0 % in microns
        offsetHeight = 0 % in microns
        
        mirroring = true; % mirror image when reaching an edge.
        
        onlineAnalysis = 'none'
        numberOfAverages = uint16(5) % number of epochs to queue
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
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % DOVES resolution  
              
            % Identify corresponding RF Filter
            [RFFilter,~,sizing] = calculateFilter(obj);
            
            % Prepare trajectory with our RF Filter.
            [traj, dist] = weightedTrajectory(obj, img2, RFFilter);
            
            % Calculate linear equivalency
            [obj.linear, obj.radii] = linearEquivalency(obj, traj, dist, sizing);
            
            % Make masks
            obj.radii = obj.radii .* (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % convert from VH to px
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            [xx, yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            m = sqrt((xx-canvasSize(1)/2).^2+(yy-canvasSize(2)/2).^2);
            
            obj.masks = zeros(canvasSize(2),canvasSize(1),obj.disks);
           
            for a = 1:size(obj.radii,2) - 1
                obj.masks(:,:,a) = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                    & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                % apply boost if desired
                if obj.boost == true
                    if obj.radii(1,a) >= obj.rig.getDevice('Stage').um2pix(obj.boostRegion(1,1))&& ...
                            obj.radii(1,a+1) <= obj.rig.getDevice('Stage').um2pix(obj.boostRegion(1,2))
                        obj.linear(a,:) = obj.linear(a,:) + obj.boostRegionBy;
                    end
                end
            end
            
            % build surround mask to ensure no pixels leak around edges.
            [xx, yy] = meshgrid(1:2*canvasSize(1),1:2*canvasSize(2));
            m = sqrt((xx-canvasSize(1)).^2+(yy-canvasSize(2)).^2);
            obj.surroundMask = m >= (max(obj.radii).*(1-obj.diskExpand/100));
            
            % Adjust axes and units for monitor (VH Pixels)
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            obj.xTraj = obj.xTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj = obj.yTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
        
            obj.counter = 1;
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('radii', obj.radii);
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
                elseif time > obj.timeTraj(end) %out of eye trajectory, hang on last frame
                    p(1) = p0(1) + obj.xTraj(end);
                    p(2) = p0(2) + obj.yTraj(end);
                else % within eye trajectory and stim time
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
            
            if obj.counter == 2
            
                % background mask (avoid leaking pixels on edge)
                surround = stage.builtin.stimuli.Rectangle();
                surround.position = canvasSize/2;
                surround.size = [2*canvasSize(1) 2*canvasSize(2)];
                surround.color = obj.backgroundIntensity;
                annulusS = uint8(obj.surroundMask*255);
                surroundMaskX = stage.core.Mask(annulusS);
                surround.setMask(surroundMaskX);
                p.addStimulus(surround);

                % disk masks
                for q = 1:size(obj.masks,3)

                    annulus = stage.builtin.stimuli.Rectangle();
                    annulus.position = canvasSize/2;
                    annulus.size = [canvasSize(1) canvasSize(2)];
                    annulusA = uint8(obj.masks(:,:,q)*255);
                    annulusMaskX = stage.core.Mask(annulusA);
                    annulus.setMask(annulusMaskX);

                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                        'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, q));

                    p.addStimulus(annulus);
                    p.addController(annulusLED);              

                end
            end
            
             % annulus background function for linear equivalent disc.
            function s = getBackground(obj, time, q)
                if time < 0 || time > obj.timeTraj(end)
                    s = obj.backgroundIntensity;
                else 
                    s = interp1(obj.timeTraj,obj.linear(q,:),time);
                end
            end
            
            obj.counter = obj.counter + 1;
            
            if obj.counter > 2
                obj.counter = 1;
            end
        end
        
        function [g, b, f] = calculateFilter(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size

            % convert to pixels
            centerSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaCenter);
            surroundSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaSurround);

            centerGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],centerSigmaPix);
            surroundGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],surroundSigmaPix);

            % calculate difference of gaussians
            diffGaussian.raw = centerGaus - surroundGaus;

            % we want to make our background light level > 0 s.t. we can
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
            
            % this parameter determines the sizing of the center/surround RF
            f(1,1) = min(centerRes);
            f(1,2) = max(centerRes);
        end
        
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
            
            % create trajectory
            for a = 1:size(obj.xTraj,2)
                % create image
               tempImg = img(obj.yTraj(1,a)-yRange:obj.yTraj(1,a)+yRange,obj.xTraj(1,a)-xRange:obj.xTraj(1,a)+xRange); 
               r(:,:,1,a) = tempImg .* tempFilt;
            end
        end
        
        function [q, radius] = linearEquivalency(obj, r, m, sizing)
            
            q = zeros(obj.disks,size(r,4));
            imgSize = [size(r,1) size(r,2)];
            
            % identify center, near-, and far-surround
            cent = 100;
            if strcmp(obj.diskFocus,'center')
                cent =  sizing(1,1) / 8;
            elseif strcmp(obj.diskFocus,'near-surround')
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
            
            maskDistribution = normpdf(1:max(imgSize)/(obj.disks):max(imgSize),cent,cent ./ obj.diskEvenness);
            maskDistribution = maskDistribution ./ (max(maskDistribution(:))) - (min(maskDistribution(:))); % normalize
            maskDistribution = abs(maskDistribution - 1);
            maskDistribution = maskDistribution ./ sum(maskDistribution(:)) .* max(imgSize);
            radius = [0 cumsum(maskDistribution)] ./ 2; % this gives us the radius of each mask
            
            if strcmp(obj.diskFocus,'none')
                radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            end
            
            for a = 1:size(radius,2) - 1
                filt = m >= radius(1,a) & m <= radius(1,a+1);
                for n = 1:size(r,4)
                    S = r(:,:,1,n) .* filt;
                    S(S==0) = [];
                    q(a,n) = mean(S) / 255;
                end
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared <= obj.numberOfAverages * 2;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted <= obj.numberOfAverages * 2;
        end
    end
end