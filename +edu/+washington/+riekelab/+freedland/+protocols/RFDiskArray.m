% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFDiskArray < edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        imageNo = 1; % natural image number (1 to 11)
        observerNo = 1; % observer number (1 to 19)
        amplification = 1; % amplify fixations by X. Setting to 0 produces a saccade-only trajectory. 
        trajectory = 'both'; % which type of stimulus to present: natural image, disk replacement, or both?        
        
        % RF field information
        rfSigmaCenter = 60; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 160; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Disk placement
        disks = 4; % number of disks that are evenly placed over a natural image.
        overrideRadii = [0 0]; % only takes effect if any value is >0. Allows any number of disks in any distribution (in px). For a 800px monitor, must contain 0 and 400 and can be represented as: [0 50 100 400].
        xSliceFrequency = 2; % how many radial slices to cut between 0 and 90 degrees.
        ySliceFrequency = 2; % how many radial slices to cut between 90 and 180 degrees.
        disksIgnoreCut = [1 2]; % starting from the center disk and moving outwards, how many disks should we NOT cut (keep circular)?

        % Disk type
        meanDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be averaged?
        contrastDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be contrasted?
        meanContrastDisks = [0 0]; % starting from the center disk and moving outwards, which disks should have both mean and contrast?
        naturalDisks = [0 0];  % starting from the center disk and moving outwards, which disks should remain a natural image?
        backgroundDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be left at background intensity?
        spatialFrequency = 100; % for contrastDisks, how many um per bar? (input from rfReverseContrastGratings) 
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        trajectoryType = symphonyui.core.PropertyType('char', 'row', {'natural','disk','both'})
        backgroundIntensity
        imageMatrix
        overrideRadiiLogical
        meanDisksLogical
        contrastDisksLogical
        xTraj
        yTraj
        timeTraj
        linear
        radii
        masks
        surroundMask
        theta
        specificOpacity
        contrast
        diskExpand
        diskOpacity
        spatialFrequencyPx
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
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',true);
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
                    
            % Identify image
            imageIdentifier = [5 8 13 17 23 27 31 39 56 64 100];
            imageVal = imageIdentifier(obj.imageNo);
            
            % Pull base trajectories and image information.
            [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal, obj.observerNo,...
                    'amplification', obj.amplification,'mirroring', true);
                
            % Scale pixels in image to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix = uint8(img2);
            
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
                
                % Set identifiers to false
                obj.meanDisksLogical = zeros(1,obj.disks);
                obj.contrastDisksLogical = zeros(1,obj.disks);
                
                % Identify corresponding RF Filter
                [RFFilter,~,~] = calculateFilter(obj);

                % Prepare trajectory with our RF Filter.
                [traj, dist] = weightedTrajectory(obj, img2, RFFilter);
                
                % Determine types of masks to calculate
                obj.meanDisks(obj.meanDisks == 0) = [];
                obj.contrastDisks(obj.contrastDisks == 0) = [];
                obj.meanContrastDisks(obj.meanContrastDisks == 0) = [];
                
                obj.meanDisksLogical(obj.meanDisks) = 1;
                obj.contrastDisksLogical(obj.contrastDisks) = 1;
                
                obj.meanDisksLogical(obj.meanContrastDisks) = 1;
                obj.contrastDisksLogical(obj.meanContrastDisks) = 1;
                
                tic

                % Calculate different disks
                if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                    [obj.linear, obj.radii, obj.contrast] = linearEquivalency(obj, traj, dist);
                else
                    [obj.linear, obj.radii, obj.theta, obj.contrast] = linearEquivalencySliced(obj, traj, dist);
                end
                
                toc

                % Recalculate masks for pixels (instead of VH pixels).
                obj.diskExpand = 0.5; % expands disks by 0.005x to ensure overlap.   
                obj.diskOpacity = 1; % opacity of disks placed.
                obj.radii = obj.radii .* (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % convert from VH to px
                
                [xx, yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
                m = sqrt((xx-canvasSize(1)/2).^2+(yy-canvasSize(2)/2).^2);
                k = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2));
                k = k ./ max(k(:)) .* 90;
                k = abs(k - 90);
                k(1:round(canvasSize(2)/2)-1,:) = k(1:round(canvasSize(2)/2)-1,:) + 180;
                
                if obj.xSliceFrequency == 0 && obj.ySliceFrequency > 0 % Special case
                    k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) = k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) + 360; % add
                end

                % The case with no slices
                if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                    
                    obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2));
                    obj.specificOpacity = ones(1,size(obj.masks,3)) .* obj.diskOpacity; % Will only place opaque masks
                    
                    for a = 1:size(obj.radii,2) - 1
                        obj.masks(:,:,a) = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                            & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                        if ismember(a,obj.naturalDisks) % If no disk...
                            obj.specificOpacity(1,a) = 0; % Make opacity zero. (Reveals natural image underneath).
                        end
                    end
                    
                else % The case with slices
                    
                    obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2),size(obj.theta,2) - 1);
                    obj.specificOpacity = ones(1,size(obj.masks,3),size(obj.masks,4)) .* obj.diskOpacity;
                    for a = 1:size(obj.radii,2) - 1
                        for b = 1:size(obj.theta,2) - 1
                            dist = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                                & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                            th = k >= (obj.theta(1,b)) & k <= (obj.theta(1,b+1));
       
                            % Ignore cuts in specific region
                            if ismember(a,obj.disksIgnoreCut)
                                th = ones(size(k));
                                obj.specificOpacity(1,a,2:end) = 0; % No need to place more than one mask.
                            end
                            
                            obj.masks(:,:,a,b) = dist .* th;
                            
                            if ismember(a,obj.naturalDisks) % If no disk...
                                obj.specificOpacity(1,a) = 0; % Make opacity zero. (Reveals natural image underneath).
                            end
                        end
                    end
                end                
            else
                obj.radii = max(canvasSize) / 2;
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
            obj.spatialFrequencyPx = obj.rig.getDevice('Stage').um2pix(obj.spatialFrequency); % Convert to pix
            if ~strcmp(obj.trajectory,'natural') 
                
                if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                    % No slices, assumes radial symmetry.
                    for q = 1:size(obj.linear,1)
                        if obj.specificOpacity(1,q) > 0 % Only place relevant disks.
                          
                            annulus = stage.builtin.stimuli.Grating();
                            annulus.position = canvasSize/2;
                            annulus.size = [canvasSize(1) canvasSize(2)]; 
                            annulus.orientation = 0; 
                            annulus.spatialFreq = 1/obj.spatialFrequencyPx;
                            annulusA = uint8(obj.masks(:,:,q)*255);
                            annulusMaskX = stage.core.Mask(annulusA);
                            annulus.setMask(annulusMaskX);

                            if ~strcmp(obj.trajectory,'both') % For single run thru.
                                if ismember(q,obj.contrastDisks) % Contrast
                                    annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                    'contrast', @(state)getContrast(obj, state.time - obj.preTime/1e3, obj.contrast(q,:)));
                                    annulus.color = obj.backgroundIntensity * 2; % No color
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLEC); % Add disk
                                    
                                elseif ismember(q,obj.meanDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, obj.linear(q,:)));
                                    annulus.contrast = 0; % No contrast
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLED); % Add disk
                                    
                                elseif ismember(q,obj.meanContrastDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, obj.linear(q,:)));
                                    annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                    'contrast', @(state)getContrast(obj, state.time - obj.preTime/1e3, obj.contrast(q,:)));
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLEC); % Add disks
                                    p.addController(annulusLED); 
                                
                                elseif ismember(q,obj.backgroundDisks)
                                    annulus.contrast = 0;
                                    annulus.color = obj.backgroundIntensity * 2;  % natively halves value for contrast resolution.
                                    p.addStimulus(annulus);
                                end
                                
                                % Only allow the scene to be visible at the exact time.
                                sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                                p.addController(sceneVisible);
                                
                            else % For successive (double length) epoch.
                                
                                cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
                                if ismember(q,obj.contrastDisks) % Contrast
                                    annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                        'contrast', @(state)getContrast(obj, state.time - cycleTime - obj.preTime/1e3, obj.contrast(q,:)));
                                    annulus.color = obj.backgroundIntensity * 2; % No color
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLEC); % Add disk
                                    
                                elseif ismember(q,obj.meanDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, obj.linear(q,:)));
                                    annulus.contrast = 0; % No contrast
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLED); % Add disk
                                    
                                elseif ismember(q,obj.meanContrastDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, obj.linear(q,:)));
                                    annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                        'contrast', @(state)getContrast(obj, state.time - cycleTime - obj.preTime/1e3, obj.contrast(q,:)));
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLEC); % Add disks
                                    p.addController(annulusLED);
                                
                                elseif ismember(q,obj.backgroundDisks)
                                    annulus.contrast = 0;
                                    annulus.color = obj.backgroundIntensity * 2; % natively halves value for contrast resolution.
                                    p.addStimulus(annulus);
                                end

                                % Only allow the scene to be visible at the exact time.
                                sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                    @(state)state.time >= cycleTime + obj.preTime * 1e-3 && state.time < cycleTime + (obj.preTime + obj.stimTime) * 1e-3);
                                p.addController(sceneVisible);
                            end              
                        end
                    end
                    
                else % With slices, more complex averaging.
                    for q = 1:size(obj.linear,1)
                        for s = 1:size(obj.linear,2)
                            if obj.specificOpacity(1,q,s) > 0
                                
                                annulus = stage.builtin.stimuli.Grating();
                                annulus.position = canvasSize/2;
                                annulus.size = [canvasSize(1) canvasSize(2)]; 
                                annulus.orientation = 0; 
                                annulus.spatialFreq = 1/obj.spatialFrequencyPx;
                                annulusA = uint8(obj.masks(:,:,q,s)*255);
                                annulusMaskX = stage.core.Mask(annulusA);
                                annulus.setMask(annulusMaskX);

                                F = obj.linear(q,s,:);
                                F = reshape(F, [1 size(F,3)]); % turn into vector

                                G = obj.contrast(q,s,:);
                                G = reshape(G, [1 size(G,3)]); % turn into vector

                                if ~strcmp(obj.trajectory,'both') % For single run thru.
                                    if ismember(q,obj.contrastDisks) % Contrast
                                        annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                        'contrast', @(state)getContrast(obj, state.time - obj.preTime/1e3, G(1,:)));
                                        annulus.color = obj.backgroundIntensity * 2; % No color
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLEC); % Add disk

                                    elseif ismember(q,obj.meanDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, F(1,:)));
                                        annulus.contrast = 0; % No contrast
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLED); % Add disk

                                    elseif ismember(q,obj.meanContrastDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, F(1,:)));
                                        annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                        'contrast', @(state)getContrast(obj, state.time - obj.preTime/1e3, G(1,:)));
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLEC); % Add disks
                                        p.addController(annulusLED); 

                                    elseif ismember(q,obj.backgroundDisks)
                                        annulus.contrast = 0;
                                        annulus.color = obj.backgroundIntensity * 2;  % natively halves value for contrast resolution.
                                        p.addStimulus(annulus);
                                    end

                                    % Only allow the scene to be visible at the exact time.
                                    sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                        @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                                    p.addController(sceneVisible);
                                
                                else % For successive (double length) epoch.
                                
                                    cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
                                    if ismember(q,obj.contrastDisks) % Contrast
                                        annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                            'contrast', @(state)getContrast(obj, state.time - cycleTime - obj.preTime/1e3, G(1,:)));
                                        annulus.color = obj.backgroundIntensity * 2; % No color
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLEC); % Add disk

                                    elseif ismember(q,obj.meanDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                            'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, F(1,:)));
                                        annulus.contrast = 0; % No contrast
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLED); % Add disk

                                    elseif ismember(q,obj.meanContrastDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                            'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, F(1,:)));
                                        annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                            'contrast', @(state)getContrast(obj, state.time - cycleTime - obj.preTime/1e3, G(1,:)));
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLEC); % Add disks
                                        p.addController(annulusLED);

                                    elseif ismember(q,obj.backgroundDisks)
                                        annulus.contrast = 0;
                                        annulus.color = obj.backgroundIntensity * 2;  % natively halves value for contrast resolution.
                                        p.addStimulus(annulus);
                                    end

                                    % Only allow the scene to be visible at the exact time.
                                    sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                        @(state)state.time >= cycleTime + obj.preTime * 1e-3 && state.time < cycleTime + (obj.preTime + obj.stimTime) * 1e-3);
                                    p.addController(sceneVisible);
                                end
                                
                            end
                        end
                    end
                end
            end
            
            % Apply far-surround mask
            surround = stage.builtin.stimuli.Rectangle();
            surround.position = canvasSize/2;
            surround.size = [2*canvasSize(1) 2*canvasSize(2)];
            surround.color = obj.backgroundIntensity;
            annulusS = uint8(obj.surroundMask*255);
            surroundMaskX = stage.core.Mask(annulusS);
            surround.setMask(surroundMaskX);
            p.addStimulus(surround);
            
            % Match mask color to linear equivalency in real time.
            function s = getBackground(obj, time, equivalency)
                if time < 0
                    s = obj.backgroundIntensity;
                elseif time > obj.timeTraj(end)
                    s = equivalency(1,end) * 2; % Hold last color
                else
                    s = interp1(obj.timeTraj,equivalency,time) * 2;
                end
            end
            
            % Match mask contrast in real time.
            function s = getContrast(obj, time, equivalency)
                if time < 0
                    s = 0;
                elseif time > obj.timeTraj(end)
                    s = equivalency(1,end); % hold last contrast
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
        function [q, radius, v] = linearEquivalency(obj, r, m)
            
            q = zeros(obj.disks,size(r,4));
            v = zeros(obj.disks,size(r,4));
            imgSize = [size(r,1) size(r,2)];
            minimumPixels = 12; % minimum amount of pixels for us to get a reasonable statistic.

            radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            
            if obj.overrideRadiiLogical == true
                radius = obj.overrideRadii ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % in VH pixels
            end
                       
            if sum([obj.meanDisksLogical,obj.contrastDisksLogical]) > 0 % Check whether the whole trajectory is needed.
                % Calculate statistics
                for a = 1:size(radius,2) - 1
                    if sum([obj.meanDisksLogical(1,a),obj.contrastDisksLogical(1,a)]) > 0 % Only regions where we absolutely need to calculate (speed increase)
                        filt = m >= radius(1,a) & m <= radius(1,a+1); % Logical mask.
                        check = filt .* m;
                        check(check == 0) = [];

                        if size(check,2) < minimumPixels % Not enough pixels to continue
                            error('Not enough pixels in region. Please change mask radii.')
                        end

                        % Apply filter across full trajectory.
                        for n = 1:size(r,4)
                            S = r(:,:,1,n) .* filt;
                            S(S==0) = [];

                            if isempty(S) % Becomes an issue when pixel weight is very small (approx 0).
                                S = 0;
                            end

                            if obj.meanDisksLogical(1,a) == 1
                                q(a,n) = mean(S) / 255;
                            end

                            if obj.contrastDisksLogical(1,a) == 1
                                v(a,n) = sqrt(var(S)) / 255;
                            end
                        end
                    end
                end
            end
        end
        
        % Calculate linear equivalency with slices
        function [q, radius, theta, v] = linearEquivalencySliced(obj, r, m)
            
            imgSize = [size(r,1) size(r,2)];
            minimumPixels = 12; % minimum amount of pixels for us to get a reasonable statistic.
            
            [xx,yy] = meshgrid(1:imgSize(2),1:imgSize(1));
            k = atan((xx - imgSize(2)/2) ./ (yy - imgSize(1)/2)); % Calculate theta
            k = k ./ max(k(:)) .* 90; % Convert to degrees
            k = abs(k - 90); % Rotate for proper polar coordinates
            k(1:floor(imgSize(1)/2),:) = k(1:floor(imgSize(1)/2),:) + 180;

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
            
            q = zeros(obj.disks,size(theta,2)-1,size(r,4)); % in form [disk #, angle, mean value]
            v = zeros(obj.disks,size(theta,2)-1,size(r,4)); % in form [disk #, angle, mean value]
            
            if sum([obj.meanDisksLogical,obj.contrastDisksLogical]) > 0 % Check whether the whole trajectory is needed.
                % Calculate statistics
                for a = 1:size(radius,2) - 1
                    for b = 1:size(theta,2) - 1
                        if sum([obj.meanDisksLogical(1,a),obj.contrastDisksLogical(1,a)]) > 0 % Only regions where we absolutely need to calculate (speed increase)
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

                            % Apply filter across full trajectory.
                            for n = 1:size(r,4)
                                S = r(:,:,1,n) .* filt;
                                S(S==0) = [];

                                if isempty(S)
                                    S = 0;
                                end

                                if obj.meanDisksLogical(1,a) == 1
                                    q(a,b,n) = mean(S) / 255;
                                end

                                if obj.contrastDisksLogical(1,a) == 1
                                    v(a,b,n) = sqrt(var(S)) / 255;
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