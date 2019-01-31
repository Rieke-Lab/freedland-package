% Amplify fixations with specific visual field.
classdef amplifyFixationsRF < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        imageNo = 1; % natural image number
        observerNo = 1; % observer number
        amplification = [0 2]; % set of fixational amplifications
        
        rfSigmaCenter = 60; % (um) Enter from difference of gaussians fit, uses to calculate center size
        rfSigmaSurround = 220; % (um) Enter from difference of gaussians fit, uses to calculate surround size
        maskLocation = 'center and surround'; % regions we place masks over
        maskStyle = 'none, lin eq disc, background'; % type of masks we place over regions.
        linearIntegrationFunction = 'gaussian'; % type of statistic for linear equivalent disc
        
        annulus = true; % solid halo between aperture and surround.
        annulusSize = 100; % beginning of surround RF, in um
        expandCenter = 0.5; % expand center outer edge by percentage to ensure mask covers effectively
        contractAnnulus = 0.5; % contract annulus inner edge by percentage to ensure mask covers effectively
        expandAnnulus = 0.5; % expand annulus outer edge by percentage to ensure mask covers effectively
        contractSurround = 0.5; % contrast surround inner edge by percentage to ensure mask covers effectively
         
        offsetWidth = 0 % in microns
        offsetHeight = 0 % in microns
        
        mirroring = true; % mirror image when reaching an edge.
        
        onlineAnalysis = 'none'
        numberOfAverages = uint16(2) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        maskLocationType = symphonyui.core.PropertyType('char', 'row', {'center', 'surround', 'center and surround'})
        maskStyleType = symphonyui.core.PropertyType('char', 'row', {'none', 'lin eq disc', 'background',...
            'none and lin eq disc', 'none and background','lin eq disc and background',...
            'none, lin eq disc, background'})
        linearIntegrationFunctionType = symphonyui.core.PropertyType('char', 'row', {'gaussian', 'uniform'})
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
            
            % Calculate center diameter
            obj.apertureDiameter = obj.rfSigmaCenter * 2 * 1.5; % we only consider 1.5sigma.
            
            % Assign variables for code efficiency.
            if obj.annulus == false
                obj.annulusSize = 0;
            end
            
            % Convert masks to pixels                
            obj.apertureDiameterPx = obj.rig.getDevice('Stage').um2pix(obj.apertureDiameter);
            obj.annulusSizePx = obj.rig.getDevice('Stage').um2pix(obj.annulusSize);
            obj.rfSigmaCenterPx = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaCenter);
            obj.rfSigmaSurroundPx = obj.rig.getDevice('Stage').um2pix(obj.rfSigmaSurround);
            
            obj.surroundRadiusPx = obj.apertureDiameterPx/2+obj.annulusSizePx;
            
            % Calculate local image statistics for lin eq disc
            if strcmp(obj.maskStyle, 'none, lin eq disc, background') ||...
                strcmp(obj.maskStyle, 'none and lin eq disc') ||...
                strcmp(obj.maskStyle, 'lin eq disc and background') || ...
                strcmp(obj.maskStyle, 'lin eq disc')
                obj.centerStat = zeros(size(obj.xTraj));
                obj.surroundStat = zeros(size(obj.xTraj));
                for n = 1:size(obj.amplification,2)
                    [cen, sur] = edu.washington.riekelab.freedland.scripts.imageSpotMean(obj.xTraj(n,:),obj.yTraj(n,:),img,obj.apertureDiameterPx,...
                        obj.annulusSizePx,obj.linearIntegrationFunction,obj.rfSigmaCenterPx,obj.rfSigmaSurroundPx);
                    obj.centerStat(n,:) = cen;
                    obj.surroundStat(n,:) = sur;
                end
            end

            % Adjust axes and units for monitor (VH Pixels)
            % We use VH pixels for our changes in motion.
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            obj.xTraj = obj.xTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj = obj.yTraj .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            
            %%%%%%%%%%%%%%%%%%%%%
            % To avoid adaptation, we randomize the order that we sample our parameters.
            % We will use a 5-digit number to code this:
            
            ampIndex = 1:size(obj.amplification,2); % amplification number in ones place
            
            % Distinguishes parameters
            if strcmp(obj.maskLocation,'center and surround')
                centerLocationIndex = 20000; %  center location to ten-thous place
                surroundLocationIndex = 200; %  center fill type to hundreds place
            elseif strcmp(obj.maskLocation,'center')
                centerLocationIndex = 20000; 
                surroundLocationIndex = 100; %  placeholder
            elseif strcmp(obj.maskLocation,'surround')
                centerLocationIndex = 10000; %  placeholder
                surroundLocationIndex = 200;
            end
            
            % Concert type of center/surround ordering
            if strcmp(obj.maskStyle, 'none, lin eq disc, background')
                maskStyleIndex = [1 2 3]; % none = 1, background = 2, lin eq disc = 3;
            elseif strcmp(obj.maskStyle, 'none and lin eq disc')
                maskStyleIndex = [1 3];
            elseif strcmp(obj.maskStyle, 'none and background')
                maskStyleIndex = [1 2];
            elseif strcmp(obj.maskStyle, 'lin eq disc and background')
                maskStyleIndex = [2 3];
            elseif strcmp(obj.maskStyle, 'none')
                maskStyleIndex = 1;
            elseif strcmp(obj.maskStyle, 'background')
                maskStyleIndex = 2;
            elseif strcmp(obj.maskStyle, 'lin eq disc')
                maskStyleIndex = 3;
            end
            
            % Code filters to only be placed in relevant locations
            if centerLocationIndex == 20000
                centerStyleIndex = maskStyleIndex * 1000; % codes to thous place
            else
                centerStyleIndex = 0;
            end
            if surroundLocationIndex == 200
                surroundStyleIndex = maskStyleIndex * 10; % codes to tens place
            else
                surroundStyleIndex = 0;
            end
            
            % We will combine each parameter to fill our space of
            % possiblities. Produces a 1xN matrix of random codes.
            A(1:size(centerLocationIndex,2),1) = centerLocationIndex;
            B(1,1:size(centerStyleIndex,2)) = centerStyleIndex;
            C(1,1,1:size(surroundLocationIndex,2)) = surroundLocationIndex;
            D(1,1,1,1:size(surroundStyleIndex,2)) =surroundStyleIndex;
            E(1,1,1,1,1:size(ampIndex,2)) = ampIndex;
            [A,B,C,D,E] = ndgrid(A,B,C,D,E); % make into repeating grid
            comb = size(centerStyleIndex,2) * size(surroundStyleIndex,2) * size(ampIndex,2);
            obj.randomIndex = reshape(A+B+C+D+E,1,comb);
            
            %%%%%%%%%%%%%%%%%%%%%
            % We will remove certain parameters from the running.
            % aperture: background, surround: background, annulus: background
            for n = 1:size(obj.randomIndex,2)
                valStr = num2str(obj.randomIndex(1,n));
                if strcmp(valStr(1:4),'2222') && obj.annulus == true
                    obj.randomIndex(1,n) = 0; % remove stimulus
                end
            end
            obj.randomIndex(obj.randomIndex == 0) = [];
            %%%%%%%%%%%%%%%%%%%%%
              
            obj.randomIndex = obj.randomIndex(randperm(size(obj.randomIndex,2))); % place in random order
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
            epoch.addParameter('microns per pixel', obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'))
                        
            obj.decoder = num2str(obj.randomIndex(obj.selectionIndex));
            fillType = {'N/A','none', 'background', 'linear equivalent disc'};

            % For analysis: if centerLocation or surroundLocation = 0, we
            % did not place a mask there.
            epoch.addParameter('centerLocationType', fillType{1,str2double(obj.decoder(2))+1});
            epoch.addParameter('surroundLocationType', fillType{1,str2double(obj.decoder(4))+1});
            epoch.addParameter('specificAmplification', obj.amplification(str2double(obj.decoder(5))));
        end
        
        function p = createPresentation(obj)
            
            % Prep stage for presentation
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            obj.distanceMatrix = createDistanceMatrix(canvasSize(1),canvasSize(2)); % make matrices of pixel distances

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
            
            % Select the random index
            ampIndex = str2double(RTdecoder(5));
            obj.xTrajSp = obj.xTraj(ampIndex,:);
            obj.yTrajSp = obj.yTraj(ampIndex,:);
            obj.centerStatSp = obj.centerStat(ampIndex,:);
            obj.surroundStatSp = obj.surroundStat(ampIndex,:);
            
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
            
            % We now have a fully formed trajectory. Let's add masks...
            
            if obj.annulus == true % Create annulus
                rect = stage.builtin.stimuli.Rectangle();
                rect.position = canvasSize/2;
                rect.size = [canvasSize(1) canvasSize(2)];
                rect.color = obj.backgroundIntensity;
                % The annulus isn't perfect, so we add 2% size as insurance.
                annulusB = uint8((obj.distanceMatrix > (obj.apertureDiameterPx/2)*(1-(obj.contractAnnulus/100)) & obj.distanceMatrix < (obj.surroundRadiusPx)*(1+(obj.expandAnnulus/100)))*255);
                centerMask = stage.core.Mask(annulusB);
                rect.setMask(centerMask);
                p.addStimulus(rect);
            end
            
            % Extract parameters from our coder         
            centerLocationIndex = str2double(RTdecoder(1));
            centerStyleIndex = str2double(RTdecoder(2));
            surroundLocationIndex = str2double(RTdecoder(3));
            surroundStyleIndex = str2double(RTdecoder(4));
            
            % Build center filter independently
            if centerLocationIndex == 2 % center
                if centerStyleIndex == 2 % background
                    disc = stage.builtin.stimuli.Ellipse();
                    disc.radiusX =  obj.apertureDiameterPx/2 * (1+(obj.expandCenter/100));
                    disc.radiusY = obj.apertureDiameterPx/2 * (1+(obj.expandCenter/100));
                    disc.color = obj.backgroundIntensity;
                    disc.position = canvasSize/2;
                    p.addStimulus(disc);
                end
                if centerStyleIndex == 3 % linear equivalent disc
                    disc = stage.builtin.stimuli.Ellipse();
                    disc.radiusX = obj.apertureDiameterPx/2 * (1+(obj.expandCenter/100));
                    disc.radiusY = obj.apertureDiameterPx/2 * (1+(obj.expandCenter/100));
                    disc.position = canvasSize/2;
                    
                    centerLED = stage.builtin.controllers.PropertyController(disc,...
                        'color', @(state)getCenterBackground(obj, state.time - obj.preTime/1e3));
                                        
                    p.addStimulus(disc);
                    p.addController(centerLED);
                end
            end
            
            %  Build surround filter independently
            if surroundLocationIndex == 2 % center
                if surroundStyleIndex == 2 % background
                    surround = stage.builtin.stimuli.Rectangle();
                    surround.position = canvasSize/2;
                    surround.color = obj.backgroundIntensity;
                    surround.size = 2.*[max(canvasSize) max(canvasSize)];
                    surroundMask = stage.core.Mask.createCircularAperture((obj.surroundRadiusPx/max(canvasSize)*(1-(obj.contractSurround/100))), 1024); %circular aperture
                    surround.setMask(surroundMask);
                    p.addStimulus(surround); %add aperture
                end
                if surroundStyleIndex == 3 % linear equivalent disc
                    surround = stage.builtin.stimuli.Rectangle();
                    surround.position = canvasSize/2;
                    surround.size = 2.*[max(canvasSize) max(canvasSize)];
                    surroundMask = stage.core.Mask.createCircularAperture((obj.surroundRadiusPx/max(canvasSize)*(1-(obj.contractSurround/100))), 1024); %circular aperture
                    surround.setMask(surroundMask);
                    
                    surroundLED = stage.builtin.controllers.PropertyController(surround,...
                        'color', @(state)getSurroundBackground(obj, state.time - obj.preTime/1e3));
                       
                    p.addStimulus(surround); %add aperture
                    p.addController(surroundLED);
                end
            end
            
            % Choose our next parameter
            if obj.selectionIndex < size(obj.randomIndex,2)
                obj.selectionIndex = obj.selectionIndex + 1; % choose next index
            elseif obj.selectionIndex == size(obj.randomIndex,2) && obj.numEpochsPrepared < obj.numberOfAverages * size(obj.randomIndex,2)
                obj.selectionIndex = 1;
                obj.randomIndex(randperm(size(obj.randomIndex,2)));
            end
            
            % Distance matrix function called for making apertures
            function m = createDistanceMatrix(sizeX,sizeY)
                [xx, yy] = meshgrid(1:sizeX,1:sizeY);
                m = sqrt((xx-(sizeX/2)).^2+(yy-(sizeY/2)).^2);
            end
            
            % center background function for linear equivalent disc.
            function c = getCenterBackground(obj, time)
                if time < 0 || time >= obj.timeTraj(end)
                    c = obj.backgroundIntensity;
                else 
                    c = interp1(obj.timeTraj,obj.centerStatSp,time);
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
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * size(obj.randomIndex,2);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * size(obj.randomIndex,2);
        end

    end
    
end