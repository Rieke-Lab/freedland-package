% Replace a natural image with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFFullFieldDiskFlash < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        imageNo = 56; % natural image number (1 to 101)
        observerNo = 1; % observer number (1 to 19) 
        frameNumber = 921; % specific frame in a eye movement trajectory. set to zero to be random.
        
        % RF field information
        rfSigmaCenter = 70; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 170; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        maskPresets = 1:12; % see draft. between 1 and 12.
        addRawMean = true; % add raw image mean as a control
        subunitIntegration = 'linear'; % currently under development. no effect.
        subunitWeight = [1 1 1 1]; % places weight by quadrant ([0 90 180 270] degrees).
        centerSurroundIntegration = 'linear'; % currently under development. no effect.
        centerSurroundWeight = [1 1 1]; % places weight on specific regions: [center annulus surround].
        
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
        disks
        rawRadius
        meanIntegration
        stat
        meanDisks
        backgroundDisks
        xSliceFrequency
        ySliceFrequency
        rotateSlices
        disksIgnoreCut
        imageMatrix
        bigCounter
        integratedStat
        surroundMask
        mainMask
        rfSizing
        meanIntensity
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
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',2); 
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));

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
        
            % Identify corresponding RF Filter
            [RFFilter,~,obj.rfSizing] = calculateFilter(obj);

            % Prepare trajectory with our RF Filter.
            [wTraj, dist, unwTraj, RFFilterVH] = weightedTrajectory(obj, img2, RFFilter);
            
            obj.meanIntensity = mean(unwTraj(:)) / 255;
            
            % Find preset properties
            obj.disks = 3;
            obj.rawRadius = [0 0.75*obj.rfSizing(1) obj.rfSizing(2) 400]; % in px
            obj.meanIntegration = 'gaussian';
            obj.stat = cell(length(obj.maskPresets),1);
            obj.integratedStat = zeros(length(obj.maskPresets),1);
            
            for a = 1:length(obj.maskPresets)
                presetNo = obj.maskPresets(a);
                if presetNo == 1        % no center mean
                    obj.meanDisks = [2 3];
                    obj.backgroundDisks = 1;
                    obj.xSliceFrequency = 0;
                    obj.ySliceFrequency = 0;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 0;
                elseif presetNo == 2    % only center mean
                    obj.meanDisks = 1;
                    obj.backgroundDisks = [2 3];
                    obj.xSliceFrequency = 0;
                    obj.ySliceFrequency = 0;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 0;
                elseif presetNo == 3    % center + annulus mean
                    obj.meanDisks = [1 2];
                    obj.backgroundDisks = 3;
                    obj.xSliceFrequency = 0;
                    obj.ySliceFrequency = 0;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 0;
                elseif presetNo == 4    % full-field mean
                    obj.meanDisks = [1 2 3];
                    obj.backgroundDisks = 0;
                    obj.xSliceFrequency = 0;
                    obj.ySliceFrequency = 0;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 0;
                elseif presetNo == 5    % center cuts only
                    obj.meanDisks = [1 2 3];
                    obj.backgroundDisks = 0;
                    obj.xSliceFrequency = 1;
                    obj.ySliceFrequency = 1;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = [2 3];
                elseif presetNo == 6    % center + annulus cuts
                    obj.meanDisks = [1 2 3];
                    obj.backgroundDisks = 0;
                    obj.xSliceFrequency = 1;
                    obj.ySliceFrequency = 1;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 3;
                elseif presetNo == 7    % full-field cuts
                    obj.meanDisks = [1 2 3];
                    obj.backgroundDisks = 0;
                    obj.xSliceFrequency = 1;
                    obj.ySliceFrequency = 1;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 0;
                elseif presetNo == 8    % center only cuts
                    obj.meanDisks = 1;
                    obj.backgroundDisks = [2 3];
                    obj.xSliceFrequency = 1;
                    obj.ySliceFrequency = 1;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = [2 3];
                elseif presetNo == 9    % center only cuts (rotated)
                    obj.meanDisks = 1;
                    obj.backgroundDisks = [2 3];
                    obj.xSliceFrequency = 1;
                    obj.ySliceFrequency = 1;
                    obj.rotateSlices = 45;
                    obj.disksIgnoreCut = [2 3];
                elseif presetNo == 10   % center only multiple cuts
                    obj.meanDisks = 1;
                    obj.backgroundDisks = [2 3];
                    obj.xSliceFrequency = 2;
                    obj.ySliceFrequency = 2;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = [2 3];
                elseif presetNo == 11   % no annulus mean
                    obj.meanDisks = [1 3];
                    obj.backgroundDisks = 2;
                    obj.xSliceFrequency = 0;
                    obj.ySliceFrequency = 0;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 0;
                elseif presetNo == 12   % no annulus cuts
                    obj.meanDisks = [1 2 3];
                    obj.backgroundDisks = 0;
                    obj.xSliceFrequency = 1;
                    obj.ySliceFrequency = 1;
                    obj.rotateSlices = 0;
                    obj.disksIgnoreCut = 2;
                end
                
                obj.stat{a,1} = findStatistic(obj, wTraj, dist, RFFilterVH, unwTraj);
                
                if ismember(presetNo,[1 2 3 4 11]) % no subunits
                    temp = obj.stat{a,1} .* obj.centerSurroundWeight; % add weight
                    obj.integratedStat(a,1) = sum(temp) / sum(obj.centerSurroundWeight);
                elseif ismember(presetNo,[5 6 7 8 9 12]) % one subunit
                    temp = obj.stat{a,1} .* repmat(obj.subunitWeight,1,3);
                    centerSurround = [nansum(temp(1:4))/sum(obj.subunitWeight)...
                        nansum(temp(5:8))/sum(obj.subunitWeight) nansum(temp(9:12))/sum(obj.subunitWeight) ];
                    obj.integratedStat(a,1) = sum(centerSurround) / sum(obj.centerSurroundWeight);
                elseif presetNo == 10 % 2 subunits
                    subunitWeighting = repelem(obj.subunitWeight,2);
                    temp = obj.stat{a,1} .* repmat(subunitWeighting,1,3);
                    centerSurround = [nansum(temp(1:8))/sum(obj.subunitWeight)...
                        nansum(temp(9:16))/sum(subunitWeighting) nansum(temp(17:24))/sum(subunitWeighting) ];
                    obj.integratedStat(a,1) = sum(centerSurround) / sum(obj.centerSurroundWeight);
                end
            end
            
            if obj.addRawMean == true
                obj.stat = [{obj.meanIntensity}; obj.stat];
                obj.integratedStat = [obj.meanIntensity; obj.integratedStat];
            end
            
            obj.imageMatrix = uint8(unwTraj);          

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
            
            obj.bigCounter = 1;
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = 2 * (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('radii', obj.rawRadius); % in pixels
            epoch.addParameter('rfSize',obj.rfSizing); % in pixels
            epoch.addParameter('frameNumber',obj.frameNumber); % in pixels
            epoch.addParameter('presetNo',obj.maskPresets(obj.bigCounter)); % in pixels
            epoch.addParameter('individualStats',obj.stat{obj.bigCounter});
            epoch.addParameter('singleMean',obj.integratedStat(obj.bigCounter)); 
            epoch.addParameter('meanDisks',obj.meanDisks);
            epoch.addParameter('backgroundDisks',obj.backgroundDisks);
            epoch.addParameter('xSliceFrequency',obj.xSliceFrequency);
            epoch.addParameter('ySliceFrequency',obj.ySliceFrequency);
            epoch.addParameter('rotateSlices', obj.rotateSlices);
            epoch.addParameter('disksIgnoreCut',obj.disksIgnoreCut);
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
            
            if obj.obj.addRawMean == true
                obj.bigCounter = mod(obj.bigCounter,length(obj.maskPresets)+1) + 1;
            else
                obj.bigCounter = mod(obj.bigCounter,length(obj.maskPresets)) + 1;
            end
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
            mask.color = obj.integratedStat(obj.bigCounter);
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
        
        % calculate value
        function val = findStatistic(obj, r, m, RFFilterVH, unwTraj)
           
            radius = obj.rawRadius ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % in VH pixels
            
            % The case with no cuts
            if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                val = zeros(1,size(radius,2) - 1);
                for a = 1:size(radius,2) - 1
                    if ismember(a,obj.meanDisks)
                        filt = m >= radius(1,a) & m <= radius(1,a+1); % Logical mask.
                        regionSize = sum(filt(:)); % pixels in area

                        if strcmp(obj.meanIntegration,'gaussian')
                            T = RFFilterVH .* filt;
                            R = sum(T(:)) / regionSize * 255; % Mean equivalence
                            S = r .* filt; % RF weighted trajectory
                        else
                            S = unwTraj .* filt; % Unweighted trajectory
                        end

                        if isempty(S) % Becomes an issue when pixel weight is very small (approx 0).
                            S = 0;
                        end

                        if strcmp(obj.meanIntegration,'gaussian')
                            val(1,a) = sum(S(:)) / (regionSize  * R); % Renormalize
                        else
                            val(1,a) = sum(S(:)) / (regionSize * 255); % Raw average
                        end
                    else
                        val(1,a) = obj.backgroundIntensity;
                    end
                end
                
            % Case with cuts
            else
                             
                imgSize = size(r);            
                [xx,yy] = meshgrid(1:imgSize(2),1:imgSize(1));
                k = atan((xx - imgSize(2)/2) ./ (yy - imgSize(1)/2)); % Calculate theta
                k = k ./ max(k(:)) .* 90; % Convert to degrees
                k = abs(k - 90); % Rotate for proper polar coordinates
                k(1:floor(imgSize(1)/2),:) = k(1:floor(imgSize(1)/2),:) + 180;
                k(imgSize(1)/2,1:imgSize(2)/2) = 180; % adjust
                k(imgSize(1)/2,imgSize(2)/2:end) = 0; % adjust
                k = mod(k - obj.rotateSlices,360); % rotate as needed.
                
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

                val = zeros(1,(size(radius,2) - 1) * (size(theta,2) - 1));
                
                counter = 1;
                
                for a = 1:size(radius,2) - 1
                    skipDisk = 0;
                    for b = 1:size(theta,2) - 1
                        if ismember(a,obj.meanDisks)
                            radiusFilt = m >= radius(1,a) & m <= radius(1,a+1); % Radial filter (r)
                            angFilt = k >= theta(1,b) & k <= theta(1,b+1); % Angular filter (theta)
                            
                            if ismember(a,obj.disksIgnoreCut) % Ignore angular filter
                                angFilt = ones(size(angFilt));
                                if b > 1
                                    skipDisk = 1;
                                end
                            end  
                            
                            filt = radiusFilt .* angFilt;
                            regionSize = sum(filt(:)); % pixels in area

                            if strcmp(obj.meanIntegration,'gaussian')
                                T = RFFilterVH .* filt;
                                R = sum(T(:)) / regionSize * 255; % Mean equivalence
                                S = r .* filt; % RF weighted trajectory
                            else
                                S = unwTraj .* filt; % Unweighted trajectory
                            end

                            if isempty(S) % Becomes an issue when pixel weight is very small (approx 0).
                                S = 0;
                            end

                            if strcmp(obj.meanIntegration,'gaussian')
                                val(1,counter) = sum(S(:)) / (regionSize  * R); % Renormalize
                            else
                                val(1,counter) = sum(S(:)) / (regionSize * 255); % Raw average
                            end
                            
                            if skipDisk == 1
                                val(1,counter) = NaN;
                            end

                        else
                            val(1,counter) = obj.backgroundIntensity;
                        end
                        
                        counter = counter + 1;
                    end
                end
            end
        end

        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.maskPresets);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.maskPresets);
        end
    end
end