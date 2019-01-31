% Amplify fixations with a carefully probed visual field.
classdef RFCancellation < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        imageNo = 1; % natural image number
        observerNo = 1; % observer number
        amplification = [0 2]; % set of fixational amplifications
        
        rfSigmaCenter = 60; % (um) Enter from difference of gaussians fit, uses to calculate center size
        rfSigmaSurround = 160; % (um) Enter from difference of gaussians fit, uses to calculate center size
        compressVertically = 1; % form ellipse by compressing up/down direction. makes horizontal ellipse.
        compressHorizontally = 1; % form ellipse by compressing left/right direction. makes vertical ellipse.
        rotation = 0; % rotation of ellipse (counterclockwise) in deg.
        thresholdSizing = 5; % threshold that defines mask radii
        centerMask = 'none, background, filter, complement, linear'; % type of mask over region
        annulusMask = 'none, background, filter, complement, linear'; % type of mask over region
        surroundMask = 'none, background, filter, complement, linear'; % type of mask over region
        maskOpacity = 1; % between 0 and 1
        
        expandCenter = 0.5; % expand center outer edge by percentage to ensure mask covers effectively
        contractSurround = 0.5; % contrast surround inner edge by percentage to ensure mask covers effectively
         
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
        centerMaskType = symphonyui.core.PropertyType('char', 'row', {'none', 'background',...
            'filter','complement','linear',...
            'none, background','none, filter','none, complement','none, linear'...
            'background, filter','background, complement','background, linear',...
            'filter, complement','filter, linear','complement, linear',...
            'none, background, filter','none, background, complement','none, background, linear',...
            'none, filter, complement','none, filter, linear','none, complement, linear',...
            'background, filter, complement','background, filter, linear','background, complement, linear',...
            'filter, complement, linear',...
            'none, background, filter, complement','none, background, filter, linear',...
            'none, background, complement, linear','none, filter, complement, linear',...
            'background, filter, complement, linear',...
            'none, background, filter, complement, linear'})
        annulusMaskType = symphonyui.core.PropertyType('char', 'row', {'none', 'background',...
            'filter','complement','linear',...
            'none, background','none, filter','none, complement','none, linear'...
            'background, filter','background, complement','background, linear',...
            'filter, complement','filter, linear','complement, linear',...
            'none, background, filter','none, background, complement','none, background, linear',...
            'none, filter, complement','none, filter, linear','none, complement, linear',...
            'background, filter, complement','background, filter, linear','background, complement, linear',...
            'filter, complement, linear',...
            'none, background, filter, complement','none, background, filter, linear',...
            'none, background, complement, linear','none, filter, complement, linear',...
            'background, filter, complement, linear',...
            'none, background, filter, complement, linear'})
        surroundMaskType = symphonyui.core.PropertyType('char', 'row', {'none', 'background',...
            'filter','complement','linear',...
            'none, background','none, filter','none, complement','none, linear'...
            'background, filter','background, complement','background, linear',...
            'filter, complement','filter, linear','complement, linear',...
            'none, background, filter','none, background, complement','none, background, linear',...
            'none, filter, complement','none, filter, linear','none, complement, linear',...
            'background, filter, complement','background, filter, linear','background, complement, linear',...
            'filter, complement, linear',...
            'none, background, filter, complement','none, background, filter, linear',...
            'none, background, complement, linear','none, filter, complement, linear',...
            'background, filter, complement, linear',...
            'none, background, filter, complement, linear'})
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON','OFF'})
        backgroundIntensity
        imageMatrix
        xTraj
        yTraj
        xTrajSp
        yTrajSp
        timeTraj
        currentStimSet
        distanceMatrix
        randomIndex
        selectionIndex
        decoder
        annulusSizePx
        apertureDiameterPx
        rfSigmaCenterPx
        rfSigmaSurroundPx
        centerStat
        surroundStat
        annulusStat
        centerStatSp
        surroundStatSp
        screenSize
        apertureDiameter
        surroundRadiusPx
        centerRadius
        annulusRadius
        surroundRadius
        filter
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
            
            % Grab image information: do this beforehand to save memory!
            [~, ~, ~,pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal, obj.observerNo,...
                 'amplification', 0, 'mirroring', obj.mirroring);
            
            % Pre-allocate memory
            baseMovementSequence.x = zeros(size(obj.amplification,2),size(pictureInformation.saccadeTracking,2));
            baseMovementSequence.y = zeros(size(obj.amplification,2),size(pictureInformation.saccadeTracking,2));
            fixMovementSequence.x = zeros(size(obj.amplification,2),size(pictureInformation.saccadeTracking,2));
            fixMovementSequence.y = zeros(size(obj.amplification,2),size(pictureInformation.saccadeTracking,2));
            
            % Convert um to pixels
            obj.offsetHeight = obj.rig.getDevice('Stage').um2pix(obj.offsetHeight);
            obj.offsetWidth = obj.rig.getDevice('Stage').um2pix(obj.offsetWidth);
            
            % Make paths for every amplification.
            for a = 1:size(obj.amplification,2)
                [~, baseMovement, fixMovement, ~] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal, obj.observerNo,...
                    'amplification', obj.amplification(1,a),'offSetHeight', obj.offsetHeight,'offSetWidth',...
                    obj.offsetWidth, 'mirroring', obj.mirroring);
                baseMovementSequence.x(a,:) = baseMovement.x;
                baseMovementSequence.y(a,:) = baseMovement.y;
                fixMovementSequence.x(a,:) = fixMovement.x;
                fixMovementSequence.y(a,:) = fixMovement.y;
            end
                
            % Adjust image
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix = uint8(img2);
            
            % Recombine trajectories
            obj.xTraj = baseMovementSequence.x + fixMovementSequence.x;
            obj.yTraj = baseMovementSequence.y + fixMovementSequence.y;
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % DOVES resolution
            
            % Split RF field into center, surround, and annulus
            [RFFilter,~] = calculateFilter(obj);
            if strcmp(obj.cellClass,'ON')
                RFFilter = 1 - RFFilter; % create complementary filter
            end
            
            % isolate filter regions
            [cen, ann, sur] = edu.washington.riekelab.freedland.scripts.imageSpotMeanExp(obj.xTraj(1,:),obj.yTraj(1,:),img,...
                    RFFilter,obj,'linearEquiv',false);
            obj.filter.annulus.background = ann.filt;
            obj.filter.surround.background = sur.safety;
            obj.filter.center.background = cen.filt;
            
            obj.filter.annulus.gradient = ann.gradient;
            obj.filter.surround.gradient = sur.safetyGradient;
            obj.filter.center.gradient = cen.gradient;

            % Adjust axes and units for monitor (VH Pixels)
            % We use VH pixels for our changes in motion.
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            obj.xTraj = obj.xTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj = obj.yTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            
            %%%%%%%%%%%%%%%%%%%%%
            % To avoid adaptation, we randomize the order that we sample our parameters.
            % We will use a 4-digit number to code this:
            ampIndex = 1:size(obj.amplification,2); % amplification
            
            fullMaskStyle = cell({obj.centerMask;obj.annulusMask;obj.surroundMask});
            fullMaskEncode = cell(3,1);
            
            % Identify mask styles.
            for a = 1:3
                if strcmp(fullMaskStyle{a,1}, 'none, background, cancellation')
                    fullMaskEncode{a,1} = [1 2 3];
                elseif strcmp(fullMaskStyle{a,1}, 'none, background')
                    fullMaskEncode{a,1} = [1 2];
                elseif strcmp(fullMaskStyle{a,1}, 'none, cancellation')
                    fullMaskEncode{a,1} = [1 3];
                elseif strcmp(fullMaskStyle{a,1}, 'background, cancellation')
                    fullMaskEncode{a,1} = [2 3];
                elseif strcmp(fullMaskStyle{a,1}, 'none')
                    fullMaskEncode{a,1} = 1;
                elseif strcmp(fullMaskStyle{a,1}, 'background')
                    fullMaskEncode{a,1} = 2;
                elseif strcmp(fullMaskStyle{a,1}, 'cancellation')
                    fullMaskEncode{a,1} = 3;
                end
            end
              
            % assign each code a 10s place.
            centerEncode = fullMaskEncode{1,1} * 10^3;
            annulusEncode = fullMaskEncode{2,1} * 10^2;
            surroundEncode = fullMaskEncode{3,1} * 10;
            
            % We will combine each parameter to fill our space of
            % possiblities. Produces a 1xN matrix of random codes.
            A(1:size(centerEncode,2),1) = centerEncode; % each independent parameter gets a dimension
            B(1,1:size(annulusEncode,2)) = annulusEncode;
            C(1,1,1:size(surroundEncode,2)) = surroundEncode;
            D(1,1,1,1:size(ampIndex,2)) = ampIndex;
            
            [A,B,C,D] = ndgrid(A,B,C,D); % make into repeating grid
            E = A + B + C + D; % add to make a base 10 encoded number
            obj.randomIndex = E(:);

            obj.randomIndex = obj.randomIndex(randperm(size(obj.randomIndex,1))); % place in random order
            obj.selectionIndex = 1; % walks along random order
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('currentStimSet', obj.currentStimSet);
                        
            obj.decoder = num2str(obj.randomIndex(obj.selectionIndex));
            fillType = {'none', 'background', 'cancellation'};

            epoch.addParameter('centerMask', fillType{1,str2double(obj.decoder(1))});
            epoch.addParameter('annulusMask', fillType{1,str2double(obj.decoder(2))});
            epoch.addParameter('surroundMask', fillType{1,str2double(obj.decoder(3))});
            epoch.addParameter('specificAmplification', obj.amplification(str2double(obj.decoder(4))));
            
            epoch.addParameter('centerRadius', obj.rig.getDevice('Stage').pix2um(obj.centerRadius.raw));
            epoch.addParameter('surroundRadius', obj.rig.getDevice('Stage').pix2um(obj.surroundRadius.raw));
            epoch.addParameter('annulusRadius', obj.rig.getDevice('Stage').pix2um(obj.annulusRadius.raw));
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
            
            % Choose a random set of parameters and pull behavior
            RTdecoder = num2str(obj.randomIndex(obj.selectionIndex));
            
            %%%%%% Apply eye trajectories
            % Select the random index
            ampIndex = str2double(RTdecoder(4));
            obj.xTrajSp = obj.xTraj(ampIndex,:);
            obj.yTrajSp = obj.yTraj(ampIndex,:);
            
            % Apply eye trajectories to move image around
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3, p0));
            
            function p = getScenePosition(obj, time, p0)
                if time < 0
                    p = p0;
                elseif time > obj.timeTraj(end) %out of eye trajectory, hang on last frame
                    p(1) = p0(1) + obj.xTrajSp(end);
                    p(2) = p0(2) + obj.yTrajSp(end);
                else % within eye trajectory and stim time
                    dx = interp1(obj.timeTraj,obj.xTrajSp,time);
                    dy = interp1(obj.timeTraj,obj.yTrajSp,time);
                    p(1) = p0(1) + dx;
                    p(2) = p0(2) + dy;
                end
            end
            
            p.addStimulus(scene);
            p.addController(scenePosition);

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            %%%%%% Apply masks for background
            % We now have a fully formed trajectory. Let's add masks for
            % background
            
            % Extract parameters from our coder         
            centerStyleIndex = str2double(RTdecoder(1));
            annulusStyleIndex = str2double(RTdecoder(2));
            surroundStyleIndex = str2double(RTdecoder(3));
            
            % center mask
            if centerStyleIndex == 2 % background
                center = stage.builtin.stimuli.Rectangle();
                center.position = canvasSize/2;
                center.size = [canvasSize(1) canvasSize(2)];
                center.color = obj.backgroundIntensity;
                annulusC = uint8(obj.filter.center.background*255);
                centerMaskX = stage.core.Mask(annulusC);
                center.setMask(centerMaskX);
                p.addStimulus(center);
                
                % turn off after stimulus completes.
                maskCVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskCVisible);
                
            elseif centerStyleIndex == 3 % composite
                center = stage.builtin.stimuli.Rectangle();
                center.position = canvasSize/2;
                center.size = [canvasSize(1) canvasSize(2)];
                center.opacity = obj.maskOpacity;
                center.color = 0;
                weightC = uint8(obj.filter.center.gradient*255);
                annulusC = stage.core.Mask(weightC);
                center.setMask(annulusC);
                p.addStimulus(center);
                
                maskCVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskCVisible);
                
            end
             
            % annulus mask
            if annulusStyleIndex == 2 % background
                annulus = stage.builtin.stimuli.Rectangle();
                annulus.position = canvasSize/2;
                annulus.size = [canvasSize(1) canvasSize(2)];
                annulus.color = obj.backgroundIntensity;
                annulusA = uint8(obj.filter.annulus.background*255);
                annulusMaskX = stage.core.Mask(annulusA);
                annulus.setMask(annulusMaskX);
                p.addStimulus(annulus);
                
                maskAVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskAVisible);
                
            elseif annulusStyleIndex == 3 % composite
                annulus = stage.builtin.stimuli.Rectangle();
                annulus.position = canvasSize/2;
                annulus.size = [canvasSize(1) canvasSize(2)];
                annulus.opacity = obj.maskOpacity;
                annulus.color = 0;
                weightA = uint8(obj.filter.annulus.gradient*255);
                annulusA = stage.core.Mask(weightA);
                annulus.setMask(annulusA);
                p.addStimulus(annulus);
                
                maskAVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskAVisible);
                
            end
            
            % surround mask
            if surroundStyleIndex == 2 % background
                surround = stage.builtin.stimuli.Rectangle();
                surround.position = canvasSize/2;
                surround.size = [2*canvasSize(1) 2*canvasSize(2)];
                surround.color = obj.backgroundIntensity;
                annulusS = uint8(obj.filter.surround.background*255);
                surroundMaskX = stage.core.Mask(annulusS);
                surround.setMask(surroundMaskX);
                p.addStimulus(surround);
                
                maskSVisible = stage.builtin.controllers.PropertyController(surround, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskSVisible);
                
            elseif surroundStyleIndex == 3 % composite
                surround = stage.builtin.stimuli.Rectangle();
                surround.position = canvasSize/2;
                surround.size = [2*canvasSize(1) 2*canvasSize(2)];
                surround.opacity = obj.maskOpacity;
                surround.color = 0;
                weightS = uint8(obj.filter.surround.gradient*255);
                annulusS = stage.core.Mask(weightS);
                surround.setMask(annulusS);
                p.addStimulus(surround);
                
                maskSVisible = stage.builtin.controllers.PropertyController(surround, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskSVisible);
                
            end
            
            % Choose our next parameter
            if obj.selectionIndex < size(obj.randomIndex,1)
                obj.selectionIndex = obj.selectionIndex + 1; % choose next index
            elseif obj.selectionIndex == size(obj.randomIndex,1) && obj.numEpochsPrepared < obj.numberOfAverages * size(obj.randomIndex,1)
                obj.selectionIndex = 1;
                obj.randomIndex(randperm(size(obj.randomIndex,1)));
            end 
        end
        
        function [g, b] = calculateFilter(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size

            % convert to pixels
            centerSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaCenter);
            surroundSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaSurround);

            GC = fspecial('gaussian',[2*canvasSize(2)*obj.compressVertically 2*canvasSize(1)*obj.compressHorizontally],centerSigmaPix);
            GR = fspecial('gaussian',[2*canvasSize(2)*obj.compressVertically 2*canvasSize(1)*obj.compressHorizontally],surroundSigmaPix);
            
            % calculate difference of gaussians
            DoG = GC - GR;
            
            % create ellipsivity
            d = imresize(DoG,[2*canvasSize(2) 2*canvasSize(1)],'bilinear');
            dR = imrotate(d,obj.rotation);
            gR = dR(round(size(dR,1)/2 - canvasSize(2)/2):round(size(dR,1)/2 + canvasSize(2)/2),...
                round(size(dR,2)/2 - canvasSize(1)/2):round(size(dR,2)/2 + canvasSize(1)/2));
            
            % we want to make our background light level > 0 s.t. we can
            % observe inhibitory effects
            DoG = gR ./ max(gR(:)); % pre-normalize filter
            a = abs(min(DoG(:))); % pre-normalize background intensity
            
            DoG = DoG+a; % scaled  s.t. no zeros exist
            
            g = DoG./max(DoG(:)); % re-normalize filter
            b = a / max(DoG(:)); % re-normalize background intensity
            
            %%%% SLIGHT DEVIATION %%%%
            % We can use this moment to quickly calculate the size of our
            % RFs.
            f = b ./ obj.thresholdSizing; % threshold by which we define our annulus.
            h = ones(1,size(g,2)) .* f;
            
            % We choose 1D and compare subtract our threshold s.t. values of 0 define the boundary.            
            res = abs(g(round(canvasSize(2)/2),:) - h);
            [~, J] = sort(res); % sort by distance
            
            % take the top few values
            K = sort(J(1,1:canvasSize(2)/10));
            placeholder = 1;
            counter = 1;
            F = zeros(4,1);
            
            for v = 2:canvasSize(2)/10
                if abs(K(1,v) - K(1,v-1)) > 10
                    F(counter,1) = median(K(1,placeholder:v-1)); % find median
                    counter = counter + 1;
                    placeholder = v;
                end
            end
            
            F(4,1) = median(K(1,placeholder:v)); % This holds four important values
            F = abs(F-canvasSize(1)/2); % center
            
            obj.centerRadius.raw = abs(mean([abs(F(2)),abs(F(3))]));
            obj.surroundRadius.raw = abs(mean([abs(F(1)),abs(F(4))]));
            obj.annulusRadius.raw = abs(mean([abs(F(1)-F(2)),abs(F(4)-F(3))]));
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * size(obj.randomIndex,1);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * size(obj.randomIndex,1);
        end
    end
    
end