% Compares full-field naturalistic image to a uniform disk
% By J. Freedland, 2019.
classdef RFFlashEquivalency < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        imageNo = 5; % natural image number (1 to 101)
        frameNumber = 600; % specific frame in a eye movement trajectory. set to zero to be random.
        observerNo = 1; % observer number (1 to 19) 
        resolution = 0.1; % luminance increase each step (0 to 1)
        minMaxIntensity = [0 0.5]; % [minimum maximum] luminance presented (0 to 1) 
        randomize = true; % randomize each luminance presented
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        subunitIntegrationType = symphonyui.core.PropertyType('char', 'row', {'linear','nonlinear'})
        centerSurroundIntegrationType = symphonyui.core.PropertyType('char', 'row', {'linear','nonlinear'})
        backgroundIntensity
        xTraj
        yTraj
        imageMatrix
        surroundMask
        mainMask
        luminance
        counter
        rfSigmaCenter
        rfSigmaSurround
        selection
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)

            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            
            % create luminance trajectory
            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.rfFlashFigureRatio',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,...
                'stimTime',obj.stimTime,'tailTime',obj.tailTime,'xval','luminance'); 

            % Pull base trajectories and image information.
            [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.imageNo, obj.observerNo,...
                    'amplification', 1,'mirroring', true);
                
            % Scale pixels in image to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            
            % Produce trajectories
            obj.xTraj = baseMovement.x + fixMovement.x;
            obj.yTraj = baseMovement.y + fixMovement.y;
            
            % Choose random frame
            if obj.frameNumber == 0
                obj.frameNumber = randperm(length(obj.xTraj));
                obj.frameNumber = obj.frameNumber(1);
            end
            obj.xTraj = obj.xTraj(1,obj.frameNumber);
            obj.yTraj = obj.yTraj(1,obj.frameNumber);

            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Calculate screen size
        
            % Overly complex way to pull individual frame (formatted from
            % other protocols).
            obj.rfSigmaCenter = 70;
            obj.rfSigmaSurround = 170;
            [RFFilter,~,obj.rfSizing] = calculateFilter(obj);
            [~, ~, unwTraj, ~] = weightedTrajectory(obj, img2, RFFilter);
            obj.imageMatrix = uint8(unwTraj); % pull individual frame    

            % There may be leaky pixels around the edge that could make
            % comparisons difficult. So, we build an mask to keep comparison controlled.
            [xx, yy] = meshgrid(1:2*canvasSize(1),1:2*canvasSize(2));
            m = sqrt((xx-canvasSize(1)).^2+(yy-canvasSize(2)).^2);
            obj.surroundMask = m >= max(canvasSize) / 2;
            protectiveMask = zeros(2*canvasSize(2),2*canvasSize(1));
            protectiveMask(round(canvasSize(2) - ceil(canvasSize(2)/2) + 1) : round(canvasSize(2) + ceil(canvasSize(2)/2) - 1), ...
                round(canvasSize(1) - ceil(canvasSize(1)/2)) + 1 : round(canvasSize(1) + ceil(canvasSize(1)/2)) - 1) = 1;
            protectiveMask = abs(protectiveMask - 1);
            obj.surroundMask = obj.surroundMask + protectiveMask; 
            
            obj.mainMask = m <= max(canvasSize);
            
            obj.luminance = obj.minMaxIntensity(1):obj.resolution:obj.minMaxIntensity(2);
            obj.counter = 1;
            
            if obj.randomize == true
                obj.selection = randperm(length(obj.luminance));
            else
                obj.selection = 1:length(obj.luminance);
            end
            
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = 2 * (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('luminance', obj.luminance(obj.selection(obj.counter)));

            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
            
            obj.counter = mod(obj.counter,length(obj.luminance)) + 1;
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();             
            p = stage.core.Presentation(2 * (obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);

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
            
            p.addStimulus(scene);
            
            % Apply eye trajectories to move image around.
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
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
            
            % Apply center mask
            mask = stage.builtin.stimuli.Rectangle();
            mask.position = canvasSize/2;
            mask.size = [canvasSize(1) canvasSize(2)];
            mask.color = obj.luminance(obj.selection(obj.counter));
            annulusM = uint8(obj.mainMask*255);
            surroundMaskM = stage.core.Mask(annulusM);
            mask.setMask(surroundMaskM);
            p.addStimulus(mask);
            
            sceneVisible2 = stage.builtin.controllers.PropertyController(mask, 'visible', ...
                @(state)state.time >= cycleTime + obj.preTime * 1e-3 && state.time < cycleTime + (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible2);
            
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
            
            % create image
            u = img(centeredYTraj-yRange:centeredYTraj+yRange-1,...
                centeredXTraj-xRange:centeredXTraj+xRange-1); 
            r = u .* tempFilt;
            
            o = tempFilt;
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.luminance);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.luminance);
        end
    end
end