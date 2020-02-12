% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFDiskArray < edu.washington.riekelab.freedland.protocols.RFDiskArrayProtocol
    
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
        disks = 3; % number of disks to replace image with.
        overrideRadii = [0 0.75 2 3]; % only takes effect if any value is >0. Arranges disks in any distribution depending on coordinate system. For info on RF coordinate system, see RFConversion function in code.
        overrideCoordinate = 'RF'; % type of coordinates to measure disk radii.
        xSliceFrequency = 1; % how many radial slices to cut between 0 and 90 degrees.
        ySliceFrequency = 1; % how many radial slices to cut between 90 and 180 degrees.
        rotateSlices = 0; % degrees to rotate all slices (0 to 90).
        disksIgnoreCut = [0 0]; % starting from the center disk and moving outwards, how many disks should we NOT cut (keep circular)?

        % Disk type
        meanDisks = [1 2 3]; % starting from the center disk and moving outwards, which disks should be averaged?
        naturalDisks = [0 0];  % starting from the center disk and moving outwards, which disks should remain a natural image?
        backgroundDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be left at background intensity?
        switchDisks = [0 0]; % starting from the center disk and moving outwards, which disks should switch intensity based on fixation?
        meanIntegration = 'gaussian'; % type of linear integration
        
        % For sequence
        playSequence = true; % plays a predefined sequence of protocols
        randomize = true; % plays sequence of protocols in random order
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(8) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        trajectoryType = symphonyui.core.PropertyType('char', 'row', {'natural','disk','both'})
        meanIntegrationType = symphonyui.core.PropertyType('char', 'row', {'uniform','gaussian'})
        overrideCoordinateType = symphonyui.core.PropertyType('char', 'row', {'pixels','RF'})
        backgroundIntensity
        imageMatrix
        overrideRadiiLogical
        meanDisksLogical
        xTraj
        yTraj
        timeTraj
        linear
        radii
        masks
        surroundMask
        theta
        specificOpacity
        diskExpandRadius
        diskExpandAngle
        diskOpacity
        rfSizing
        preUnit
        postUnit
        switchTraj
        switchVal
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.freedland.protocols.RFDiskArrayProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)

            prepareRun@edu.washington.riekelab.freedland.protocols.RFDiskArrayProtocol(obj);
            
            % Must have prerender set to avoid lagging through the replacement trajectory.     
            if ~strcmp(obj.trajectory,'natural') && obj.rig.getDevice('Stage').getConfigurationSetting('prerender') == 0
                error('Must have prerender set') 
            end
            
            % Define settings as needed.
            redefineSettings(obj, true);

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

            % Catch common errors.
            checkDisks = sum(1:obj.disks);
            checkAssignments = sum([obj.meanDisks obj.switchDisks obj.naturalDisks obj.backgroundDisks]);
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
                
            % Scale image pixels to monitor
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix = uint8(img2);

            % Produce trajectories
            obj.xTraj = baseMovement.x + fixMovement.x;
            obj.yTraj = baseMovement.y + fixMovement.y;
            
            % Establish switching intensities
            if sum(obj.switchDisks) > 0
                saccades = find(pictureInformation.saccadeTracking == 1);
                switchLocations = double(diff(saccades) > 1) .* saccades(2:end);
                switchLocations(switchLocations == 0) = [];
                switchLocations = [1 switchLocations length(obj.xTraj)];
                
                if obj.switchVal == 1
                    intensities = [0 obj.backgroundIntensity.*2]; % [min, max] intensity between disks
                elseif obj.switchVal == 2
                    intensities = [obj.backgroundIntensity.*2 0]; % swap order
                end
                
                obj.switchTraj = zeros(1,length(obj.xTraj));
                
                counter = 0;
                for a = 1:length(switchLocations)-1
                    obj.switchTraj(switchLocations(a):switchLocations(a+1)) = intensities(counter+1);
                    counter = mod(counter+1,2);
                end
            end
            
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
                [RFFilter,obj.rfSizing] = calculateFilter(obj);

                % Prepare trajectory with our RF Filter
                [wTraj, unwTraj, RFFilterVH] = weightedTrajectory(obj, img2, RFFilter);

                % For an override, convert RF units into pixels.
                if strcmp(obj.overrideCoordinate,'RF') && obj.overrideRadiiLogical == true                   
                    RFConversion(obj)
                end

                % Calculate different disk statistics (all in VH units)
                tic
                if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                    [obj.linear, obj.radii] = linearEquivalency(obj, wTraj, unwTraj, RFFilterVH);
                else
                    [obj.linear, obj.radii, obj.theta] = linearEquivalencySliced(obj, wTraj, unwTraj, RFFilterVH);
                end
                toc

                %%% We have calculated all statistics in DOVES (VH) units. Now, we need to
                %%% define the space (in pixels) that masks will be placed over for our monitor.
                
                % Recalculate masks for pixels (instead of VH pixels).
                obj.diskExpandRadius = 0.5;   % expands disks radially by percentage to ensure overlap.  
                obj.diskExpandAngle = 0; % expands disks' angle coverage by percentage to ensure overlap.  
                obj.diskOpacity = 1;    % opacity of disks placed. This is a good setting to play with during stimulus testing. (Matches pixels underneath).
                
                obj.preUnit = obj.radii;
                changeUnits(obj,'VH2PIX'); % convert from VH (DOVES) to px
                obj.radii = obj.postUnit;
                
                % In pixels
                [xx, yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
                m = sqrt((xx-canvasSize(1)/2).^2+(yy-canvasSize(2)/2).^2);
                k = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2));
                k = k ./ max(k(:)) .* 90;
                k = abs(k - 90);
                k(1:round(canvasSize(2)/2)-1,:) = k(1:round(canvasSize(2)/2)-1,:) + 180;
                
                k = mod(k - obj.rotateSlices,360);
                
                if obj.xSliceFrequency == 0 && obj.ySliceFrequency > 0 % Special case
                    k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) = k(round(canvasSize(2)/2):end,round(canvasSize(1)/2):end) + 360; % add
                end

                % The case with no slices
                if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                    
                    obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2));
                    obj.specificOpacity = ones(1,size(obj.masks,3)) .* obj.diskOpacity; % Will only place opaque masks
                    
                    for a = 1:size(obj.radii,2) - 1
                        obj.masks(:,:,a) = m >= (obj.radii(1,a).*(1-obj.diskExpandRadius/100))...
                            & m <= (obj.radii(1,a+1).*(1+obj.diskExpandRadius/100));
                        if ismember(a,obj.naturalDisks) % If no disk...
                            obj.specificOpacity(1,a) = 0; % Make opacity zero. (Reveals natural image underneath).
                        end
                    end
                    
                else % The case with slices
                    
                    obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2),size(obj.theta,2) - 1);
                    obj.specificOpacity = ones(1,size(obj.masks,3),size(obj.masks,4)) .* obj.diskOpacity;
                    for a = 1:size(obj.radii,2) - 1
                        for b = 1:size(obj.theta,2) - 1
                            dist = m >= (obj.radii(1,a).*(1-obj.diskExpandRadius/100))...
                                & m <= (obj.radii(1,a+1).*(1+obj.diskExpandRadius/100));
                            th = k >= (obj.theta(1,b)) & k <= (obj.theta(1,b+1));
       
                            % Add in disk expanding for radial coordinates
                            expansionCoef = 360 * (obj.diskExpandAngle/100);
                            expansionDisk1 = k >= (obj.theta(1,b) - expansionCoef) & k <= (obj.theta(1,b));
                            expansionDisk2 = k <= (obj.theta(1,b+1) + expansionCoef) & k >= (obj.theta(1,b+1));
                            
                            % Special case
                            if obj.theta(1,b) == 0
                                expansionDisk1 = k >= (360 - expansionCoef) & k <= 360;
                            end
                            if obj.theta(1,b+1) == 360
                                expansionDisk2 = k <= (0 + expansionCoef) & k >= 0;
                            end
                                
                            % Add expanded angles to disk
                            th = (th + expansionDisk1 + expansionDisk2) > 0;
                            
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
            obj.surroundMask = m >= (max(obj.radii).*(1-obj.diskExpandRadius/100));
            protectiveMask = zeros(2*canvasSize(2),2*canvasSize(1));
            protectiveMask(round(canvasSize(2) - ceil(canvasSize(2)/2) + 1) : round(canvasSize(2) + ceil(canvasSize(2)/2) - 1), ...
                round(canvasSize(1) - ceil(canvasSize(1)/2)) + 1 : round(canvasSize(1) + ceil(canvasSize(1)/2)) - 1) = 1;
            protectiveMask = abs(protectiveMask - 1);
            obj.surroundMask = obj.surroundMask + protectiveMask; 
            
            % Adjust axes and units for monitor (DOVES VH units to pixels)
            obj.xTraj = -(obj.xTraj - size(img,2)/2);
            obj.yTraj = (obj.yTraj - size(img,1)/2);
            
            obj.preUnit = obj.xTraj;
            changeUnits(obj,'VH2PIX'); % convert from VH to px
            obj.xTraj = obj.postUnit;

            obj.preUnit = obj.yTraj;
            changeUnits(obj,'VH2PIX'); % convert from VH to px
            obj.yTraj = obj.postUnit;
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.freedland.protocols.RFDiskArrayProtocol(obj, epoch);
            
            if obj.numEpochsCompleted > 0 && mod(obj.numEpochsCompleted,obj.numberOfAverages) == 0
                redefineSettings(obj, false)
            end
            
            device = obj.rig.getDevice(obj.amp);
            if ~strcmp(obj.trajectory,'both')
                duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            else
                duration = 2 * (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            end
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('radii', obj.radii); % in pixels
            epoch.addParameter('rfSize',obj.rfSizing); % in pixels
            
            if obj.playSequence == true % add metadata 
                epoch.addParameter('specificImageNo', obj.imageNo);
                epoch.addParameter('specificXSliceFrequency', obj.xSliceFrequency);
                epoch.addParameter('specificYSliceFrequency', obj.ySliceFrequency);
                epoch.addParameter('specificRotateSlices', obj.rotateSlices);
                epoch.addParameter('specificDisksIgnoreCut', obj.disksIgnoreCut);
                epoch.addParameter('specificMeanDisks', obj.meanDisks);
                epoch.addParameter('specificBackgroundDisks', obj.backgroundDisks);
                epoch.addParameter('specificNaturalDisks', obj.naturalDisks);
                epoch.addParameter('specificSwitchDisks', obj.switchDisks);
                epoch.addParameter('switchDisksMarker', obj.switchVal);
            else
                epoch.addParameter('backgroundIntensity', obj.backgroundIntensity); % not accurate for sequence
            end
            
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

            if sum(obj.naturalDisks) > 0
                % Repeat naturalistic stimulus underneath masks
                sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                    @(state)rem(state.time,cycleTime) >= obj.preTime * 1e-3 && rem(state.time,cycleTime) < (obj.preTime + obj.stimTime) * 1e-3);
            else
                sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            end
            p.addController(sceneVisible);

            %%%%%% Apply masks %%%%%% 
            if ~strcmp(obj.trajectory,'natural') 
                
                if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                    % No slices, assumes radial symmetry.
                    for q = 1:size(obj.linear,1)
                        if obj.specificOpacity(1,q) > 0 % Only place relevant disks.
                          
                            annulus = stage.builtin.stimuli.Grating();
                            annulus.position = canvasSize/2;
                            annulus.size = [canvasSize(1) canvasSize(2)]; 
                            annulus.orientation = 0; 
                            annulus.spatialFreq = 0;
                            annulus.opacity = obj.specificOpacity(1,q);
                            annulusA = uint8(obj.masks(:,:,q)*255);
                            annulusMaskX = stage.core.Mask(annulusA);
                            annulus.setMask(annulusMaskX);

                            if ~strcmp(obj.trajectory,'both') % For single run thru.
                                    
                                if ismember(q,obj.meanDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, obj.linear(q,:)));
                                    annulus.contrast = 0; % No contrast
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLED); % Add disk
                                
                                elseif ismember(q,obj.backgroundDisks)
                                    annulus.contrast = 0;
                                    annulus.color = obj.backgroundIntensity * 2;  % Natively halves value for contrast resolution.
                                    p.addStimulus(annulus);
                                    
                                elseif ismember(q,obj.switchDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, obj.switchTraj));
                                    annulus.contrast = 0; % No contrast
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLED); % Add disk
                                    
                                end
                                
                                % Only allow the scene to be visible at the exact time.
                                sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                                p.addController(sceneVisible);
                                
                            else % For successive (double length) epoch.
                                
                                cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
                                    
                                if ismember(q,obj.meanDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, obj.linear(q,:)));
                                    annulus.contrast = 0; % No contrast
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLED); % Add disk  
                                
                                elseif ismember(q,obj.backgroundDisks)
                                    annulus.contrast = 0;
                                    annulus.color = obj.backgroundIntensity * 2; % natively halves value for contrast resolution.
                                    p.addStimulus(annulus);
                                    
                                elseif ismember(q,obj.switchDisks)
                                    annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, obj.switchTraj));
                                    annulus.contrast = 0; % No contrast
                                    p.addStimulus(annulus); % Add stimulus
                                    p.addController(annulusLED); % Add disk
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
                                annulus.spatialFreq = 0;
                                annulus.opacity = obj.specificOpacity(1,q,s);
                                annulusA = uint8(obj.masks(:,:,q,s)*255);
                                annulusMaskX = stage.core.Mask(annulusA);
                                annulus.setMask(annulusMaskX);

                                F = obj.linear(q,s,:);
                                F = reshape(F, [1 size(F,3)]); % turn into vector

                                if ~strcmp(obj.trajectory,'both') % For single run thru.

                                    if ismember(q,obj.meanDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, F(1,:)));
                                        annulus.contrast = 0; % No contrast
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLED)

                                    elseif ismember(q,obj.backgroundDisks)
                                        annulus.contrast = 0;
                                        annulus.color = obj.backgroundIntensity * 2;  % natively halves value for contrast resolution.
                                        p.addStimulus(annulus);
                                        
                                    elseif ismember(q,obj.switchDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - obj.preTime/1e3, obj.switchTraj));
                                        annulus.contrast = 0; % No contrast
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLED); % Add disk    
                                    end

                                    % Only allow the scene to be visible at the exact time.
                                    sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                        @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                                    p.addController(sceneVisible);
                                
                                else % For successive (double length) epoch.
                                
                                    cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;

                                    if ismember(q,obj.meanDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                            'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, F(1,:)));
                                        annulus.contrast = 0; % No contrast
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLED); % Add disk

                                    elseif ismember(q,obj.backgroundDisks)
                                        annulus.contrast = 0;
                                        annulus.color = obj.backgroundIntensity * 2;  % natively halves value for contrast resolution.
                                        p.addStimulus(annulus);
                                        
                                    elseif ismember(q,obj.switchDisks)
                                        annulusLED = stage.builtin.controllers.PropertyController(annulus,...
                                        'color', @(state)getBackground(obj, state.time - cycleTime - obj.preTime/1e3, obj.switchTraj));
                                        annulus.contrast = 0; % No contrast
                                        p.addStimulus(annulus); % Add stimulus
                                        p.addController(annulusLED); % Add disk
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
                    s = obj.backgroundIntensity .* 2;
                elseif time > obj.timeTraj(end)
                    s = equivalency(1,end) * 2; % Hold last color
                else
                    s = interp1(obj.timeTraj,equivalency,time) * 2;
                end
            end
        end
        
        % Calculate RF Filter
        function [normFilter, f] = calculateFilter(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Identify screen size

            % Convert to pixels
            obj.preUnit = obj.rfSigmaCenter;
            changeUnits(obj,'UM2PIX');
            centerSigmaPix = obj.postUnit;
            
            obj.preUnit = obj.rfSigmaSurround;
            changeUnits(obj,'UM2PIX');
            surroundSigmaPix = obj.postUnit;

            % Generate 2D gaussians
            centerGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],centerSigmaPix);
            surroundGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],surroundSigmaPix);

            % Calculate difference of gaussians
            diffGaussian = centerGaus - surroundGaus;
            normFilter = diffGaussian ./ max(diffGaussian(:)); % Normalize filter

            % Shift gaussian s.t. all values are positive
            shift = abs(min(normFilter(:)));  
            normFilter = normFilter + shift;
            normFilter = normFilter ./ max(normFilter(:)); % Renormalize

            % Take 2D slice along longest axis.
            if size(normFilter,1) > size(normFilter,2)
                slice = normFilter(round(median(1:size(normFilter,1)):size(normFilter,1)),round(median(1:size(normFilter,2))));
            else
                slice = normFilter(round(median(1:size(normFilter,1))),round(median(1:size(normFilter,2)):size(normFilter,2)));
            end
            threshold = 0.02; % Tends to work nicely for many RFs

            sizing = find(slice <= threshold); % Find annulus
            
            % This parameter determines the sizing of the center/surround RF
            f(1,1) = min(sizing); % in pixels
            f(1,2) = max(sizing);
        end
        
        % Apply RF Filter over the entire trajectory.
        function [w, u, m] = weightedTrajectory(obj, img, RFFilter)

            % Calculate trajectory size
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            
            obj.preUnit = canvasSize;
            changeUnits(obj,'PIX2VH'); % VH units for all DOVES movies.
            imgSize = floor(obj.postUnit);
            
            % We want an even amount of pixels on each side of the trajectory.
            xLength = floor(imgSize(1) / 2);
            yLength = floor(imgSize(2) / 2);
            
            % Calculate pixel range for every frame.
            xRange = zeros(length(obj.xTraj),xLength.*2+1);
            yRange = zeros(length(obj.yTraj),yLength.*2+1);
            for a = 1:length(obj.xTraj)
                xRange(a,:) = obj.xTraj(a) - xLength : obj.xTraj(a) + xLength;
                yRange(a,:) = obj.yTraj(a) - yLength : obj.yTraj(a) + yLength;
            end
            
            % In Stage, we will enlarge our image to fit the monitor.
            % For our calculations, we will do the opposite. By "shrinking" 
            % our monitor to fit the image, we improve computation time considerably.
            m = imresize(RFFilter, [size(yRange,2) size(xRange,2)]);
            
            % Create trajectory  
            w = zeros(size(yRange,2),size(xRange,2),1,length(obj.xTraj));
            u = zeros(size(yRange,2),size(xRange,2),1,length(obj.xTraj));
            for a = 1:length(obj.xTraj)
               % Pull image
               u(:,:,1,a) = img(yRange(a,:),xRange(a,:)); 
               
               % Convolve image with RF
               w(:,:,1,a) = u(:,:,1,a) .* m;
            end
        end
        
        % Calculate linear equivalency without slices
        function [q, radius] = linearEquivalency(obj, wTraj, unwTraj, RFFilterVH)
            
            q = zeros(obj.disks,size(wTraj,4));
            imgSize = [size(wTraj,1) size(wTraj,2)];
            
            [xx, yy] = meshgrid(1:imgSize(2),1:imgSize(1));
            m = sqrt((xx-imgSize(2)/2).^2+(yy-imgSize(1)/2).^2);
            
            minimumPixels = 8; % Minimum amount of pixels for us to get a reasonable statistic.

            radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2; % Evenly distributed (in VH)
            
            if obj.overrideRadiiLogical == true
                obj.preUnit = obj.overrideRadii;
                changeUnits(obj,'PIX2VH'); % in VH pixels
                radius = obj.postUnit;
            end
                       
            if sum(obj.meanDisksLogical) > 0 % Check whether the whole trajectory is needed.
                
                % Calculate statistics
                for a = 1:size(radius,2) - 1
                    if obj.meanDisksLogical(1,a) > 0 % Only regions where we absolutely need to calculate (speed increase)
                        
                        filt = m >= radius(1,a) & m <= radius(1,a+1); % Logical mask (VH coordinates)
                        regionSize = sum(filt(:)); % Total number of pixels
                        
                        if regionSize < minimumPixels % Not enough pixels to continue
                            error('Not enough pixels in region. Please change mask radii.')
                        end

                        if strcmp(obj.meanIntegration,'gaussian')
                            T = RFFilterVH .* filt;             % Individual region's RF
                            R = sum(T(:)) / regionSize * 255;   % Normalizing factor for region
                        end

                        % Apply filter across full trajectory.
                        for n = 1:size(wTraj,4)
                            
                            if obj.meanDisksLogical(1,a) == 1 % Disk belongs
                                
                                if strcmp(obj.meanIntegration,'gaussian')
                                    S = wTraj(:,:,1,n) .* filt;     % RF weighted trajectory
                                else
                                    S = unwTraj(:,:,1,n) .* filt;   % Unweighted trajectory
                                end
                                
                                if isempty(S) % Becomes an issue when pixel weight is very small (approx 0).
                                    S = 0;
                                end
                                
                                if strcmp(obj.meanIntegration,'gaussian')
                                    q(a,n) = sum(S(:)) / (regionSize  * R);     % Renormalize
                                else
                                    q(a,n) = sum(S(:)) / (regionSize * 255);    % Raw average
                                end
                                
                            end

                        end
                    end
                end
            end
        end
        
        % Calculate linear equivalency with slices
        function [q, radius, theta] = linearEquivalencySliced(obj, wTraj, unwTraj, RFFilterVH)
            
            imgSize = [size(wTraj,1) size(wTraj,2)];
            minimumPixels = 8; % minimum amount of pixels for us to get a reasonable statistic.
            
            % Define regions for sliced disks
            [xx,yy] = meshgrid(1:imgSize(2),1:imgSize(1));
            m = sqrt((xx-imgSize(2)/2).^2+(yy-imgSize(1)/2).^2);
            
            k = atan((xx - imgSize(2)/2) ./ (yy - imgSize(1)/2)); % Calculate theta
            k = k ./ max(k(:)) .* 90;                             % Convert to degrees
            k = abs(k - 90);                                      % Rotate for proper polar coordinates
            k(1:floor(imgSize(1)/2),:) = k(1:floor(imgSize(1)/2),:) + 180;
            
            k = mod(k - obj.rotateSlices,360); % Rotate accordingly

            radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            
            if obj.overrideRadiiLogical == true
                obj.preUnit = obj.overrideRadii;
                changeUnits(obj,'PIX2VH'); % in VH pixels
                radius = obj.postUnit;
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
            
            q = zeros(obj.disks,size(theta,2)-1,size(wTraj,4)); % in form [disk #, angle, mean value]
            
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
                            regionSize = sum(filt(:));

                            if regionSize < minimumPixels % Not enough pixels to continue
                                error('Not enough pixels in region. Please change mask radii.')
                            end

                            if strcmp(obj.meanIntegration,'gaussian')
                                T = RFFilterVH .* filt;
                                R = sum(T(:)) / regionSize * 255; % Mean equivalence
                            end

                            % Apply filter across full trajectory.
                            for n = 1:size(wTraj,4)

                                if obj.meanDisksLogical(1,a) == 1 % Disk belongs

                                    if strcmp(obj.meanIntegration,'gaussian')
                                        S = wTraj(:,:,1,n) .* filt; % RF weighted trajectory
                                    else
                                        S = unwTraj(:,:,1,n) .* filt; % Unweighted trajectory
                                    end

                                    if isempty(S) % Becomes an issue when pixel weight is very small (approx 0).
                                        S = 0;
                                    end

                                    if strcmp(obj.meanIntegration,'gaussian')
                                        q(a,b,n) = sum(S(:)) / (regionSize  * R); % Renormalize
                                    else
                                        q(a,b,n) = sum(S(:)) / (regionSize  * 255); % Raw average
                                    end

                                end

                            end
                        end
                    end
                end
            end
        end
        
        function changeUnits(obj,type)
            
            if strcmp(type,'UM2PIX')
                % um / (um/pix) = pix
                obj.postUnit = obj.preUnit ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
                
            elseif strcmp(type,'PIX2UM')
                % pix * (um/pix) = um
                obj.postUnit = obj.preUnit .* obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
   
            elseif strcmp(type,'UM2VH')
                % From DOVES database: 1 VH pixel = 1 arcmin = 3.3 um on monkey retina
                % um / (3.3 um/VH) = VH
                obj.postUnit = obj.preUnit ./ 3.3;
            
            elseif strcmp(type,'VH2UM')
                % VH * (3.3 um/VH) = um
                obj.postUnit = obj.preUnit .* 3.3;
            
            elseif strcmp(type,'PIX2VH')
                % (3.3 um/VH) / (um/pix) = pix/VH
                ratio = 3.3 ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
                
                % pix / (pix/VH) = VH
                obj.postUnit = obj.preUnit ./ ratio;

            elseif strcmp(type,'VH2PIX')
                % (3.3 um/VH) / (um/pix) = pix/VH
                ratio = 3.3 ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
                
                % VH * (pix/VH) = pix
                obj.postUnit = obj.preUnit .* ratio;
            else
                error('incorrect unit conversion.')
            end
        end
        
        function RFConversion(obj)
            
            % RF coordinates are a special form of units.
            % It uses a cell's RF to define borders between the
            % center, near surround (annulus), and far surround.
            % Then, it uses these borders to draw disks.
            %
            % RF coordinate values indicate the following:
            %   1 == border between center and annulus
            %   2 == border between annulus and surround
            %
            % An RF value of 0.75 or 1.75 would draw disks 25% inward from
            % each region's border.
            %
            % This coordinate system is used because it is invariant of
            % each cell's individual RF. In other words, it scales
            % accordingly for cells with large RFs and small RFs.
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Identify screen size
            
            centerAnnulusBorder = obj.rfSizing(1);
            annulusSurroundBorder = obj.rfSizing(2);
            
            centerSize = centerAnnulusBorder;
            annulusSize = annulusSurroundBorder - centerAnnulusBorder;
            surroundSize = (max(canvasSize) / 2) - annulusSurroundBorder;
            
            for a = 1:obj.disks+1
                if obj.overrideRadii(1,a) <= 1
                    RFCoordinate = obj.overrideRadii(1,a);
                    obj.overrideRadii(1,a) = RFCoordinate .* centerSize;
                elseif obj.overrideRadii(1,a) <= 2
                    RFCoordinate = obj.overrideRadii(1,a) - 1;
                    obj.overrideRadii(1,a) = centerSize + (RFCoordinate .* annulusSize);
                elseif obj.overrideRadii(1,a) <= 3
                    RFCoordinate = obj.overrideRadii(1,a) - 2;
                    obj.overrideRadii(1,a) = centerSize + annulusSize + (RFCoordinate .* surroundSize);
                end
            end
        end
        
        function redefineSettings(obj,t)
            if obj.playSequence == true
                redefineSettings@edu.washington.riekelab.freedland.protocols.RFDiskArrayProtocol(obj,t);
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            if obj.playSequence == true
                tf = obj.numEpochsPrepared < obj.numberOfAverages*obj.numberOfExperiments;
            else
                tf = obj.numEpochsPrepared < obj.numberOfAverages;
            end  
        end
        
        function tf = shouldContinueRun(obj)
            if obj.playSequence == true
                tf = obj.numEpochsCompleted < obj.numberOfAverages*obj.numberOfExperiments;
            else
                tf = obj.numEpochsCompleted < obj.numberOfAverages;
            end
        end
    end
end