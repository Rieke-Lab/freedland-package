% Amplify fixations with a carefully probed visual field.
classdef RFNaturalTrajectory < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        imageNo = 1; % natural image number (1 to 11)
        observerNo = 1; % observer number (1 to 19)
        amplification = [0 2]; % set of fixational amplifications
        
        rfSigmaCenter = 60; % (um) enter from difference of gaussians fit, uses to calculate center size
        rfSigmaSurround = 160; % (um) enter from difference of gaussians fit, uses to calculate center size
        compressVertically = 1; % form ellipse by compressing up/down direction. makes horizontal ellipse.
        compressHorizontally = 1; % form ellipse by compressing left/right direction. makes vertical ellipse.
        rotation = 0; % rotation of ellipse (counterclockwise) in deg.
        centerMask = 'none, background, complement, linear'; % type of mask over region
        annulusMask = 'none, background, linear'; % type of mask over region
        surroundMask = 'none, background, linear'; % type of mask over region
        thresholdSizing = 5; % threshold that defines mask radii
        
        contractLinEqCenter = 0; % samples pixels for linear equivalent disc from % smaller area. prevents skewing from low-value pixels in annulus.
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
        annulusStatSp
        screenSize
        apertureDiameter
        surroundRadiusPx
        centerRadius
        surroundRadius
        filter
        linearMarker
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
                      
            %%%%%%%%%%%%%%%%%%%%% 
            % We can now calculate the filters for our trajectory
            
            % Split RF field into center, surround, and annulus
            [RFFilter,~] = calculateFilter(obj);
            RFFilterC = 1 - RFFilter; % complementary filter
            fullMaskStyle = cell({obj.centerMask;obj.annulusMask;obj.surroundMask});
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            
            [cen, ann, sur] = edu.washington.riekelab.freedland.scripts.imageSpotMean(obj.xTraj(1,:),obj.yTraj(1,:),img,canvasSize,...
                    RFFilter,obj,'linearEquiv',false,'umPerPix',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
                
            % background filter
            obj.filter.annulus.background = ann.filt;
            obj.filter.surround.background = sur.safety;
            obj.filter.center.background = cen.filt;
            
            % rf filter
            obj.filter.annulus.filter = ann.gradient;
            obj.filter.surround.filter = sur.safetyGradient;
            obj.filter.center.filter = cen.gradient;

            % complementary filter
            if ~isempty(strfind(fullMaskStyle{1,1},'complement')) || ~isempty(strfind(fullMaskStyle{2,1},'complement')) || ...
                    ~isempty(strfind(fullMaskStyle{3,1},'complement'))
                [cen, ann, sur] = edu.washington.riekelab.freedland.scripts.imageSpotMean(obj.xTraj(1,:),obj.yTraj(1,:),img,canvasSize,...
                        RFFilterC,obj,'linearEquiv',false); % we use RFFilterC instead of RFFilter
                obj.filter.annulus.complement = ann.gradient;
                obj.filter.surround.complement = sur.safetyGradient;
                obj.filter.center.complement = cen.gradient;
            end
            
            % linear equivalent
            if ~isempty(strfind(fullMaskStyle{1,1},'linear')) || ~isempty(strfind(fullMaskStyle{2,1},'linear')) || ...
                    ~isempty(strfind(fullMaskStyle{3,1},'linear'))
                obj.linearMarker = true;
                obj.filter.center.linear = zeros(size(obj.xTraj));
                obj.filter.surround.linear = zeros(size(obj.xTraj));
                obj.filter.annulus.linear = zeros(size(obj.xTraj));
                for n = 1:size(obj.amplification,2)
                    [cen, sur, ann] = edu.washington.riekelab.freedland.scripts.imageSpotMean(obj.xTraj(n,:),obj.yTraj(n,:),img,canvasSize,...
                        RFFilter,obj,'linearEquiv',true);
                    obj.filter.center.linear(n,:) = cen.data;
                    obj.filter.surround.linear(n,:) = sur.data;
                    obj.filter.annulus.linear(n,:) = ann.data;
                end
            else
                obj.linearMarker = false;
            end
            %%%%%%%%%%%%%%%%%%%%%

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
            fullMaskEncode = cell(3,1);
            
            % Identify mask styles.
            for a = 1:3
                if ~isempty(strfind(fullMaskStyle{a,1},'none'))
                    fullMaskEncode{a,1} = [fullMaskEncode{a,1} 1];
                end
                if ~isempty(strfind(fullMaskStyle{a,1},'background'))
                    fullMaskEncode{a,1} = [fullMaskEncode{a,1} 2];
                end
                if ~isempty(strfind(fullMaskStyle{a,1},'filter'))
                    fullMaskEncode{a,1} = [fullMaskEncode{a,1} 3];
                end
                if ~isempty(strfind(fullMaskStyle{a,1},'complement'))
                    fullMaskEncode{a,1} = [fullMaskEncode{a,1} 4];
                end
                if ~isempty(strfind(fullMaskStyle{a,1},'linear'))
                    fullMaskEncode{a,1} = [fullMaskEncode{a,1} 5];
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
            %%%%%%%%%%%%%%%%%%%%%

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
                        
            obj.decoder = num2str(obj.randomIndex(obj.selectionIndex));
            fillType = {'none', 'background', 'filter', 'complement', 'linear'};

            epoch.addParameter('centerMask', fillType{1,str2double(obj.decoder(1))});
            epoch.addParameter('annulusMask', fillType{1,str2double(obj.decoder(2))});
            epoch.addParameter('surroundMask', fillType{1,str2double(obj.decoder(3))});
            epoch.addParameter('specificAmplification', obj.amplification(str2double(obj.decoder(4))));
            
            epoch.addParameter('centerRadius', obj.rig.getDevice('Stage').pix2um(obj.centerRadius.maskR));
            epoch.addParameter('surroundRadius', obj.rig.getDevice('Stage').pix2um(obj.surroundRadius.maskR));
            epoch.addParameter('linearCenterRadius', obj.rig.getDevice('Stage').pix2um(obj.centerRadius.linR));
            
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
            
            % Choose a random set of parameters and pull behavior
            RTdecoder = num2str(obj.randomIndex(obj.selectionIndex));
            
            %%%%%% Apply eye trajectories
            % Select the random index
            ampIndex = str2double(RTdecoder(4));
            obj.xTrajSp = obj.xTraj(ampIndex,:);
            obj.yTrajSp = obj.yTraj(ampIndex,:);
            
            if obj.linearMarker == true
                obj.centerStatSp = obj.filter.center.linear(ampIndex,:);
                obj.surroundStatSp = obj.filter.surround.linear(ampIndex,:);
                obj.annulusStatSp = obj.filter.annulus.linear(ampIndex,:);
            end
            
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
                
            elseif centerStyleIndex == 3 % filter
                center = stage.builtin.stimuli.Rectangle();
                center.position = canvasSize/2;
                center.size = [canvasSize(1) canvasSize(2)];
                center.opacity = 1;
                center.color = obj.backgroundIntensity;
                weightC = uint8(obj.filter.center.complement*255); % we switch! 0s in symphony = no mask.
                annulusC = stage.core.Mask(weightC);
                center.setMask(annulusC);
                p.addStimulus(center);
                
                maskCVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskCVisible);
                
            elseif centerStyleIndex == 4 % complement
                center = stage.builtin.stimuli.Rectangle();
                center.position = canvasSize/2;
                center.size = [canvasSize(1) canvasSize(2)];
                center.opacity = 1;
                center.color = obj.backgroundIntensity;
                weightC = uint8(obj.filter.center.filter*255); % we switch! 0s in symphony = no mask.
                annulusC = stage.core.Mask(weightC);
                center.setMask(annulusC);
                p.addStimulus(center);
                
                maskCVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskCVisible);
                
            elseif centerStyleIndex == 5 % linear
                center = stage.builtin.stimuli.Rectangle();
                center.position = canvasSize/2;
                center.size = [canvasSize(1) canvasSize(2)];
                center.opacity = 1;
                annulusC = uint8(obj.filter.center.background*255);
                centerMaskX = stage.core.Mask(annulusC);
                center.setMask(centerMaskX);
                    
                centerLED = stage.builtin.controllers.PropertyController(center,...
                    'color', @(state)getCenterBackground(obj, state.time - obj.preTime/1e3));
                                        
                p.addStimulus(center);
                p.addController(centerLED);
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
                
            elseif annulusStyleIndex == 3 % filter
                annulus = stage.builtin.stimuli.Rectangle();
                annulus.position = canvasSize/2;
                annulus.size = [canvasSize(1) canvasSize(2)];
                annulus.opacity = 1;
                annulus.color = obj.backgroundIntensity;
                weightA = uint8(obj.filter.annulus.complement*255);
                annulusA = stage.core.Mask(weightA);
                annulus.setMask(annulusA);
                p.addStimulus(annulus);
                
                maskAVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskAVisible);
                
            elseif annulusStyleIndex == 4 % complement
                annulus = stage.builtin.stimuli.Rectangle();
                annulus.position = canvasSize/2;
                annulus.size = [canvasSize(1) canvasSize(2)];
                annulus.opacity = 1;
                annulus.color = obj.backgroundIntensity;
                weightA = uint8(obj.filter.annulus.filter*255);
                annulusA = stage.core.Mask(weightA);
                annulus.setMask(annulusA);
                p.addStimulus(annulus);
                
                maskAVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskAVisible);
            
            elseif annulusStyleIndex == 5 % linear
                annulus = stage.builtin.stimuli.Rectangle();
                annulus.position = canvasSize/2;
                annulus.size = [canvasSize(1) canvasSize(2)];
                annulus.opacity = 1;
                annulusA = uint8(obj.filter.annulus.background*255);
                annulusMaskX = stage.core.Mask(annulusA);
                annulus.setMask(annulusMaskX);
                
                annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                    'color', @(state)getAnnulusBackground(obj, state.time - obj.preTime/1e3));
                                        
                p.addStimulus(annulus);
                p.addController(annulusLED);              
                
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
                
            elseif surroundStyleIndex == 3 % filter
                surround = stage.builtin.stimuli.Rectangle();
                surround.position = canvasSize/2;
                surround.size = [2*canvasSize(1) 2*canvasSize(2)];
                surround.opacity = 1;
                surround.color = obj.backgroundIntensity;
                weightS = uint8(obj.filter.surround.complement*255);
                annulusS = stage.core.Mask(weightS);
                surround.setMask(annulusS);
                p.addStimulus(surround);
                
                maskSVisible = stage.builtin.controllers.PropertyController(surround, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskSVisible);
                
            elseif surroundStyleIndex == 4 % complement
                surround = stage.builtin.stimuli.Rectangle();
                surround.position = canvasSize/2;
                surround.size = [2*canvasSize(1) 2*canvasSize(2)];
                surround.opacity = 1;
                surround.color = obj.backgroundIntensity;
                weightS = uint8(obj.filter.surround.filter*255);
                annulusS = stage.core.Mask(weightS);
                surround.setMask(annulusS);
                p.addStimulus(surround);
                
                maskSVisible = stage.builtin.controllers.PropertyController(surround, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(maskSVisible);
                
            elseif surroundStyleIndex == 5 % linear
                surround = stage.builtin.stimuli.Rectangle();
                surround.position = canvasSize/2;
                surround.size = [2*canvasSize(1) 2*canvasSize(2)];
                annulusS = uint8(obj.filter.surround.background*255);
                surround.opacity = 1;
                surroundMaskX = stage.core.Mask(annulusS);
                surround.setMask(surroundMaskX);
                
                surroundLED = stage.builtin.controllers.PropertyController(surround,...
                    'color', @(state)getSurroundBackground(obj, state.time - obj.preTime/1e3));
                                        
                p.addStimulus(surround);
                p.addController(surroundLED);
            end
            
            % Choose our next parameter
            if obj.selectionIndex < size(obj.randomIndex,1)
                obj.selectionIndex = obj.selectionIndex + 1; % choose next index
            elseif obj.selectionIndex == size(obj.randomIndex,1) && obj.numEpochsPrepared < obj.numberOfAverages * size(obj.randomIndex,1)
                obj.selectionIndex = 1;
                obj.randomIndex(randperm(size(obj.randomIndex,1)));
            end 
            
            % center background function for linear equivalent disc.
            function c = getCenterBackground(obj, time)
                if time < 0 || time >= obj.timeTraj(end)
                    c = obj.backgroundIntensity;
                else 
                    c = interp1(obj.timeTraj,obj.centerStatSp,time);
                end
            end
            
             % annulus background function for linear equivalent disc.
            function s = getAnnulusBackground(obj, time)
                if time < 0 || time > obj.timeTraj(end)
                    s = obj.backgroundIntensity;
                else 
                    s = interp1(obj.timeTraj,obj.annulusStatSp,time);
                end
            end
            
            % Surround background function for linear equivalent disc.
            function s = getSurroundBackground(obj, time)
                if time < 0 || time > obj.timeTraj(end)
                    s = obj.backgroundIntensity;
                else 
                    s = interp1(obj.timeTraj,obj.surroundStatSp,time);
                end
            end
        end
        
        function [g, b] = calculateFilter(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size

            % convert to pixels
            centerSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaCenter);
            surroundSigmaPix = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaSurround);

            centerGaus = fspecial('gaussian',[2*canvasSize(2)*obj.compressVertically 2*canvasSize(1)*obj.compressHorizontally],centerSigmaPix);
            surroundGaus = fspecial('gaussian',[2*canvasSize(2)*obj.compressVertically 2*canvasSize(1)*obj.compressHorizontally],surroundSigmaPix);

            % calculate difference of gaussians
            diffGaussian.raw = centerGaus - surroundGaus;

            % create ellipsivity
            diffGaussian.exp = imresize(diffGaussian.raw,[2*canvasSize(2) 2*canvasSize(1)],'bilinear');
            diffGaussian.rot = imrotate(diffGaussian.exp,obj.rotation);
            diffGaussian.size = diffGaussian.rot(round(size(diffGaussian.rot,1)/2 - canvasSize(2)/2):round(size(diffGaussian.rot,1)/2 + canvasSize(2)/2),...
            round(size(diffGaussian.rot,2)/2 - canvasSize(1)/2):round(size(diffGaussian.rot,2)/2 + canvasSize(1)/2));

            % we want to make our background light level > 0 s.t. we can
            % observe inhibitory effects
            filte = diffGaussian.size ./ max(diffGaussian.size(:)); % pre-normalize filter
            supp = abs(min(filte(:))); % pre-normalize suppressive region

            filte = filte+supp; % scaled  s.t. no zeros exist

            g = filte./max(filte(:)); % re-normalize filter
            b = supp / max(filte(:)); % re-normalize background intensity

            %%%% SLIGHT DEVIATION %%%%
            % We can use this moment to calculate the size of our RFs.
            threshold = b ./ obj.thresholdSizing; % threshold by which we define our surround.
            undoneFilter = imrotate(g,-obj.rotation); % rerotate filter to normalize to x,y plane

            sliceX = undoneFilter(round(median(1:size(undoneFilter,1))),round(median(1:size(undoneFilter,2))):end); % take a 2D slice (X)
            sliceY = undoneFilter(round(median(1:size(undoneFilter,1))):end,round(median(1:size(undoneFilter,2)))); % take a 2D slice (Y)
            sliceX(sliceX == 0) = []; % remove edges from rotation
            sliceY(sliceY == 0) = [];

            sizingX = sliceX <= threshold; % find annulus
            sizingY = sliceY <= threshold; % find annulus

            lenX = 1:size(sizingX,2);
            lenY = 1:size(sizingY,1);
            centerResX = sizingX .* lenX; % make values index
            centerResY = sizingY .* lenY';
            centerResX(centerResX == 0) = [];
            centerResY(centerResY == 0) = [];
            
            % this parameter determines the sizing of the center/surround RF
            obj.centerRadius.trueX = min(centerResX);
            obj.centerRadius.trueY = min(centerResY);
            obj.surroundRadius.trueX = max(centerResX);
            obj.surroundRadius.trueY = max(centerResY);
            
            % this parameter determines the area of the center that we
            % average pixels over. Otherwise, we will be heavily skewed by
            % a large number of low-intensity pixels at the edge.
            obj.centerRadius.linEqX = min(centerResY) * (1 - obj.contractLinEqCenter/100); 
            obj.centerRadius.linEqY = min(centerResX) * (1 - obj.contractLinEqCenter/100);
            
            % calculate radii
            obj.centerRadius.maskR = sqrt(obj.centerRadius.trueX.^2 + obj.centerRadius.trueY.^2);
            obj.surroundRadius.maskR = sqrt(obj.surroundRadius.trueX.^2 + obj.surroundRadius.trueY.^2);
            obj.centerRadius.linR = sqrt(obj.centerRadius.linEqX.^2 + obj.centerRadius.linEqY.^2);
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared <= obj.numberOfAverages * size(obj.randomIndex,1);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted <= obj.numberOfAverages * size(obj.randomIndex,1);
        end
    end
    
end