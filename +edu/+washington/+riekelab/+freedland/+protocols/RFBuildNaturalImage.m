% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFBuildNaturalImage < edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        referenceImageNo = 1; % natural image number (1 to 11) that we calculate statistics for.
        testImageNo = 2; % natural image number (1 to 11) that we apply statistics to.
        referenceObserverNo = 1; % observer number (1 to 19) that we calculate statistics for.
        testObserverNo = 1; % observer number (1 to 19) that we apply statistics to.
        amplification = 1; % amplify fixations by X. Setting to 0 produces a saccade-only trajectory. 
        totalOpacity = 0.2; % 0 to 1
        trajectory = 'controlled'; % which type of stimulus to present: compared (img1 + filtered img2) or controlled (img1 + img2 + filtered img2)   
        
        % RF field information
        rfSigmaCenter = 70; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 170; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Disk placement
        disks = 3; % number of disks that are evenly placed over a natural image.
        overrideRadii = [0 0.75 2 3]; % only takes effect if any value is >0. Allows any number of disks in any distribution. In pixels, for a 800px monitor, must contain 0 and 400 and can be represented as: [0 50 100 400]. In RF coordinates, must contain 0 and 3, where 1, 2 are the radius of the center, surround respectively. For a [center surround] = [70 170], [0 0.5 1 1.5 2 3] in RF = [0 35 70 120 170 400] in px.
        overrideCoordinate = 'RF'; % type of coordinates to measure disk radii.
        xSliceFrequency = 1; % how many radial slices to cut between 0 and 90 degrees.
        ySliceFrequency = 1; % how many radial slices to cut between 90 and 180 degrees.
        disksIgnoreCut = [3 0]; % starting from the center disk and moving outwards, how many disks should we NOT cut (keep circular)?

        % Disk type
        meanDisks = [1 2 3]; % starting from the center disk and moving outwards, which disks should be averaged?
        naturalDisks = [0 0];  % starting from the center disk and moving outwards, which disks should remain a natural image?
        backgroundDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be left at background intensity?
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(10) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        trajectoryType = symphonyui.core.PropertyType('char', 'row', {'compared','controlled'})
        meanIntegrationType = symphonyui.core.PropertyType('char', 'row', {'uniform','gaussian'})
        contrastIntegrationType = symphonyui.core.PropertyType('char', 'row', {'uniform','gaussian','contrast'})
        overrideCoordinateType = symphonyui.core.PropertyType('char', 'row', {'pixels','RF'})
        backgroundIntensity1
        backgroundIntensity2
        imageMatrix1
        imageMatrix2
        xTraj1
        yTraj1
        xTraj2
        yTraj2
        timeTraj
        ratio
        radii
        theta
        overrideRadiiLogical
        meanDisksLogical
        masks
        surroundMask
        diskOpacity
        temporalOpacity
        diskExpand
        spatialFrequencyPx
        rfSizing
        linear
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)

            prepareRun@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj);

            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            if strcmp(obj.trajectory,'compared') % Splits the epoch into two for online analysis.
                obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',2);
            else % triple epoch
                obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                    obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',3);
            end
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Catch common errors        
            if ~strcmp(obj.trajectory,'natural') && obj.rig.getDevice('Stage').getConfigurationSetting('prerender') == 0
                error('Must have prerender set') % Pre-render required, else trajectory lags and isn't accurate
            end
            
            checkDisks = sum(1:obj.disks);
            checkAssignments = sum([obj.meanDisks obj.naturalDisks obj.backgroundDisks]);
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
            
            % Nomenclature: 1 = reference, 2 = test
                    
            % Identify image
            imageIdentifier = [5 8 13 17 23 27 31 39 56 64 100];
            imageVal1 = imageIdentifier(obj.referenceImageNo);
            imageVal2 = imageIdentifier(obj.testImageNo);
            
            % Pull base trajectories and image information.
            [~, baseMovement1, fixMovement1, pictureInformation1] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal1, obj.referenceObserverNo,...
                    'amplification', obj.amplification,'mirroring', true);
                
            [~, baseMovement2, fixMovement2, pictureInformation2] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal2, obj.testObserverNo,...
                    'amplification', obj.amplification,'mirroring', true);
                
            % Scale pixels in image to monitor
            img = pictureInformation1.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity1 = mean(img(:));
            img1 = img.*255;
            obj.imageMatrix1 = uint8(img1);
            
            img = pictureInformation2.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity2 = mean(img(:));
            img2 = img.*255;
            obj.imageMatrix2 = uint8(img2);
            
            % Produce trajectories
            obj.xTraj1 = baseMovement1.x + fixMovement1.x;
            obj.yTraj1 = baseMovement1.y + fixMovement1.y;
            obj.xTraj2 = baseMovement2.x + fixMovement2.x;
            obj.yTraj2 = baseMovement2.y + fixMovement2.y;
            
            % We do not need to consider the entire trajectory, however.
            frames = round((obj.stimTime + 50) / 1000 * 200); % max # of frames in DOVES database, with 50ms cushion
            if size(obj.xTraj1,2) <= frames
                frames = size(obj.xTraj1,2);
            end
            
            % This enables a faster computation time for shorter stimTimes
            obj.xTraj1 = obj.xTraj1(1,1:frames);
            obj.yTraj1 = obj.yTraj1(1,1:frames);
            obj.xTraj2 = obj.xTraj2(1,1:frames);
            obj.yTraj2 = obj.yTraj2(1,1:frames);
            obj.timeTraj = (0:(length(obj.xTraj1)-1)) ./ 200; % DOVES resolution

            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Calculate screen size
            
            % Identify type of disks
            obj.meanDisksLogical = zeros(1,obj.disks);
            obj.meanDisks(obj.meanDisks == 0) = [];
            obj.meanDisksLogical(obj.meanDisks) = 1;

            % Identify corresponding RF Filter
            [RFFilter,~,obj.rfSizing] = calculateFilter(obj); % consistent among both

            % Prepare trajectory with our RF Filter.
            [dist, unwTraj1, ~, ~] = calculateTrajectory(obj, img1, RFFilter);

            % For an override, convert RF units into pixels.
            if strcmp(obj.overrideCoordinate,'RF') && obj.overrideRadiiLogical == true                   
                H = obj.rfSizing(2) - obj.rfSizing(1);
                G = max(canvasSize) / 2 - obj.rfSizing(2);
                for a = 1:obj.disks+1
                    if obj.overrideRadii(1,a) <= 1
                        obj.overrideRadii(1,a) = obj.overrideRadii(1,a) * obj.rfSizing(1);
                    elseif obj.overrideRadii(1,a) <= 2
                        A = obj.overrideRadii(1,a) - 1; % normalize
                        obj.overrideRadii(1,a) = obj.rfSizing(1) + A * H;
                    elseif obj.overrideRadii(1,a) <= 3
                        A = obj.overrideRadii(1,a) - 2; % normalize
                        obj.overrideRadii(1,a) = obj.rfSizing(2) + A * G;
                    end
                end
            end

            tic

            % Calculate difference between disks
            if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                [obj.ratio, obj.radii] = buildOpacity(obj, unwTraj1, dist);
            else
                [obj.ratio, obj.radii, obj.theta] = buildOpacitySlices(obj, unwTraj1, dist);
            end

            toc

            % Recalculate masks for pixels (instead of VH pixels).
            obj.diskExpand = 0; 
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

                for a = 1:size(obj.radii,2) - 1
                    obj.masks(:,:,a) = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                        & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                end

            else % The case with slices

                obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2),size(obj.theta,2) - 1);
                for a = 1:size(obj.radii,2) - 1
                    for b = 1:size(obj.theta,2) - 1
                        dist = m >= (obj.radii(1,a).*(1-obj.diskExpand/100))...
                            & m <= (obj.radii(1,a+1).*(1+obj.diskExpand/100));
                        th = k >= (obj.theta(1,b)) & k <= (obj.theta(1,b+1));

                        % Ignore cuts in specific region
                        if ismember(a,obj.disksIgnoreCut)
                            th = ones(size(k));
                        end

                        obj.masks(:,:,a,b) = dist .* th;
                    end
                end
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
            obj.xTraj1 = -(obj.xTraj1 - size(img,2)/2);
            obj.yTraj1 = (obj.yTraj1 - size(img,1)/2);
            obj.xTraj1 = obj.xTraj1 .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj1 = obj.yTraj1 .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            
            obj.xTraj2 = -(obj.xTraj2 - size(img,2)/2);
            obj.yTraj2 = (obj.yTraj2 - size(img,1)/2);
            obj.xTraj2 = obj.xTraj2 .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.yTraj2 = obj.yTraj2 .* 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'); 
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.freedland.protocols.RepeatPrerenderStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            if strcmp(obj.trajectory,'controlled')
                duration = 3 * (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            else
                duration = 2 * (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            end
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity1);
            epoch.addParameter('radii', obj.radii); % in pixels
            epoch.addParameter('rfSize',obj.rfSizing); % in pixels
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
                                    
            if strcmp(obj.trajectory,'controlled')
                p = stage.core.Presentation(3 * (obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            else % Present both stimuli in succession
                p = stage.core.Presentation(2 * (obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            end

            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity1);
            
            % Prep to display image 1
            scene1 = stage.builtin.stimuli.Image(obj.imageMatrix1);
            scene1.size = [size(obj.imageMatrix1,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.imageMatrix1,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            p0 = canvasSize/2;
            scene1.position = p0;
            
            % Use linear interpolation when scaling the image
            scene1.setMinFunction(GL.LINEAR);
            scene1.setMagFunction(GL.LINEAR);
            
            % Apply eye trajectories to move image around.
            scene1Position = stage.builtin.controllers.PropertyController(scene1,...
                'position', @(state)getScene1Position(obj, state.time - obj.preTime/1e3, p0));
            
            function p = getScene1Position(obj, time, p0)
                if time <= 0
                    p = p0;
                elseif time > obj.timeTraj(end) % Beyond eye trajectory, hang on last frame
                    p(1) = p0(1) + obj.xTraj1(end);
                    p(2) = p0(2) + obj.yTraj1(end);
                else % Within eye trajectory and stim time
                    dx = interp1(obj.timeTraj,obj.xTraj1,time);
                    dy = interp1(obj.timeTraj,obj.yTraj1,time);
                    p(1) = p0(1) + dx;
                    p(2) = p0(2) + dy;
                end
            end
            
            p.addStimulus(scene1);
            p.addController(scene1Position);

            scene1Visible = stage.builtin.controllers.PropertyController(scene1, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(scene1Visible);
            
            scene2 = stage.builtin.stimuli.Image(obj.imageMatrix2);
            scene2.size = [size(obj.imageMatrix2,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.imageMatrix2,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            scene2.position = p0;
            
            % Use linear interpolation when scaling the image
            scene2.setMinFunction(GL.LINEAR);
            scene2.setMagFunction(GL.LINEAR);
            
            % Apply eye trajectories to move image around.
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            scene2Position = stage.builtin.controllers.PropertyController(scene2,...
                'position', @(state)getScene2Position(obj, rem(state.time - obj.preTime/1e3 - cycleTime,cycleTime), p0));
            
            function p = getScene2Position(obj, time, p0)
                if time <= 0
                    p = p0;
                elseif time > obj.timeTraj(end) % Beyond eye trajectory, hang on last frame
                    p(1) = p0(1) + obj.xTraj1(end);
                    p(2) = p0(2) + obj.yTraj1(end);
                else % Within eye trajectory and stim time
                    dx = interp1(obj.timeTraj,obj.xTraj1,time);
                    dy = interp1(obj.timeTraj,obj.yTraj1,time);
                    p(1) = p0(1) + dx;
                    p(2) = p0(2) + dy;
                end
            end
            
            p.addStimulus(scene2);
            p.addController(scene2Position);

            scene2Visible = stage.builtin.controllers.PropertyController(scene2, 'visible', ...
                @(state)rem(state.time - cycleTime,cycleTime) >= obj.preTime * 1e-3 && rem(state.time - cycleTime,cycleTime) < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(scene2Visible);

            %%%%%% Apply masks %%%%%% 
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0
                % No slices, assumes radial symmetry.
                for q = 1:size(obj.ratio,1)

                    annulus = stage.builtin.stimuli.Grating();
                    annulus.position = canvasSize/2;
                    annulus.size = [canvasSize(1) canvasSize(2)]; 
                    annulusA = uint8(obj.masks(:,:,q)*255);
                    annulusMaskX = stage.core.Mask(annulusA);
                    annulus.setMask(annulusMaskX);
                    annulus.contrast = 0; 
                    annulus.opacity = obj.totalOpacity;

                    if strcmp(obj.trajectory,'compared') % For single run thru.

                        if ismember(q,obj.meanDisks)
                            annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                            'color', @(state)getColor(obj, state.time - obj.preTime/1e3 - cycleTime, obj.ratio(q,:)));
                            p.addStimulus(annulus); % Add stimulus
                            p.addController(annulusLEC); % Add disk

                        elseif ismember(q,obj.backgroundDisks)
                            annulus.color = obj.backgroundIntensity1 * 2;  % natively halves value for contrast resolution.
                            p.addStimulus(annulus);
                        end

                        % Only allow the scene to be visible at the exact time.
                        sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                            @(state)state.time >= cycleTime + obj.preTime * 1e-3 && state.time < (cycleTime + obj.preTime + obj.stimTime) * 1e-3);
                        p.addController(sceneVisible);

                    else % For triple epoch
                        
                        if ismember(q,obj.meanDisks)
                            annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                            'color', @(state)getColor(obj, state.time - obj.preTime/1e3 - 2*cycleTime, obj.ratio(q,:)));
                            annulus.contrast = 0; % No contrast
                            p.addStimulus(annulus); % Add stimulus
                            p.addController(annulusLEC); % Add disk

                        elseif ismember(q,obj.backgroundDisks)
                            annulus.color = obj.backgroundIntensity * 2; % natively halves value for contrast resolution.
                            p.addStimulus(annulus);
                        end

                        %Only allow the scene to be visible at the exact time.
                        sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                            @(state)state.time >= 2*cycleTime + obj.preTime * 1e-3 && state.time < 2*cycleTime + (obj.preTime + obj.stimTime) * 1e-3);
                        p.addController(sceneVisible);
                    end              
                end

            else % With slices, more complex averaging.
                for q = 1:size(obj.ratio,1)
                    for s = 1:size(obj.ratio,2)

                        annulus = stage.builtin.stimuli.Grating();
                        annulus.position = canvasSize/2;
                        annulus.size = [canvasSize(1) canvasSize(2)]; 
                        annulusA = uint8(obj.masks(:,:,q,s)*255);
                        annulusMaskX = stage.core.Mask(annulusA);
                        annulus.setMask(annulusMaskX);
                        annulus.contrast = 0; 
                        annulus.opacity = obj.totalOpacity;

                        F = obj.ratio(q,s,:);
                        F = reshape(F, [1 size(F,3)]); % turn into vector

                        if strcmp(obj.trajectory,'compared')

                            if ismember(q,obj.meanDisks)
                                annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                'color', @(state)getColor(obj, state.time - obj.preTime/1e3 - cycleTime, F(1,:)));
                                p.addStimulus(annulus); % Add stimulus
                                p.addController(annulusLEC); % Add disk

                            elseif ismember(q,obj.backgroundDisks)
                                annulus.contrast = 0;
                                annulus.color = obj.backgroundIntensity1 * 2;  % natively halves value for contrast resolution.
                                p.addStimulus(annulus);
                            end

                            %Only allow the scene to be visible at the exact time.
                            sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                @(state)state.time >= cycleTime + obj.preTime * 1e-3 && state.time < (cycleTime + obj.preTime + obj.stimTime) * 1e-3);
                            p.addController(sceneVisible);

                        else % For triple length epoch

                            if ismember(q,obj.meanDisks)
                                annulusLEC = stage.builtin.controllers.PropertyController(annulus,...
                                    'color', @(state)getColor(obj, state.time - 2*cycleTime - obj.preTime/1e3, F(1,:)));
                                p.addStimulus(annulus); % Add stimulus
                                p.addController(annulusLEC); % Add disk

                            elseif ismember(q,obj.backgroundDisks)
                                annulus.contrast = 0;
                                annulus.color = obj.backgroundIntensity1 * 2;  % natively halves value for contrast resolution.
                                p.addStimulus(annulus);
                            end

                            %Only allow the scene to be visible at the exact time.
                            sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                                @(state)state.time >= 2*cycleTime + obj.preTime * 1e-3 && state.time < 2*cycleTime + (obj.preTime + obj.stimTime) * 1e-3);
                            p.addController(sceneVisible);
                        end
                    end
                end
            end
            
            % Apply far-surround mask
            surround = stage.builtin.stimuli.Rectangle();
            surround.position = canvasSize/2;
            surround.size = [2*canvasSize(1) 2*canvasSize(2)];
            surround.color = obj.backgroundIntensity1;
            annulusS = uint8(obj.surroundMask*255);
            surroundMaskX = stage.core.Mask(annulusS);
            surround.setMask(surroundMaskX);
            p.addStimulus(surround);
            
            % Match mask color to linear equivalency in real time.
            function s = getColor(obj, time, equivalency)
                if time < 0
                    s = 0;
                elseif time > obj.timeTraj(end)
                    s = equivalency(1,end) * 2; % Hold last opacity
                else
                    s = interp1(obj.timeTraj,equivalency,time) * 2;
                end
            end
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
        function [m, u, o, u2] = calculateTrajectory(obj, img, RFFilter)
            
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
            u = zeros(size(m,1),size(m,2),1,size(obj.xTraj1,2));
            u2 = zeros(size(m,1),size(m,2),1,size(obj.xTraj2,2));
            for a = 1:size(obj.xTraj1,2)
                % create image
               u(:,:,1,a) = img(obj.yTraj1(1,a)-yRange:obj.yTraj1(1,a)+yRange,obj.xTraj1(1,a)-xRange:obj.xTraj1(1,a)+xRange); 
               u2(:,:,1,a) = img(obj.yTraj2(1,a)-yRange:obj.yTraj2(1,a)+yRange,obj.xTraj2(1,a)-xRange:obj.xTraj2(1,a)+xRange); 
            end
            
            o = tempFilt;
        end
        
        % Calculate linear equivalency without slices
        function [q, radius] = buildOpacity(obj, traj1, traj2, m)
            
            q = zeros(obj.disks,size(traj1,4));
            imgSize = [size(traj1,1) size(traj2,2)];
            minimumPixels = 12; % minimum amount of pixels for us to get a reasonable statistic.

            radius = [0 cumsum(repelem(max(imgSize)/(obj.disks),1,(obj.disks)))] ./ 2;
            
            if obj.overrideRadiiLogical == true
                radius = obj.overrideRadii ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')); % in VH pixels
            end
                       
            if sum(obj.meanDisksLogical) > 0 % Check whether the whole trajectory is needed.
                % Calculate statistics
                for a = 1:size(radius,2) - 1
                    if obj.meanDisksLogical(1,a) > 0 % Only regions where we absolutely need to calculate (speed increase)
                        filt = m >= radius(1,a) & m <= radius(1,a+1); % Logical mask.
                        check = filt;
                        check(check == 0) = [];

                        if size(check,2) < minimumPixels % Not enough pixels to continue
                            error('Not enough pixels in region. Please change mask radii.')
                        end

                        % Apply filter across full trajectory.
                        for n = 1:size(traj1,4)
                            
                            if obj.meanDisksLogical(1,a) == 1 % Disk belongs
                                
                                S1 = traj1(:,:,1,n) .* filt; % reference img
                                S1(S1 == 0) = [];

                                q(a,n) = mean(S1(:)) / 255;
                            end
                        end
                    end
                end
            end
        end
        
        % Calculate linear equivalency with slices
        function [q, radius, theta] = buildOpacitySlices(obj, traj1, m)
             
            imgSize = [size(traj1,1) size(traj1,2)];
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
            
            q = zeros(obj.disks,size(theta,2)-1,size(traj1,4)); % in form [disk #, angle, mean value]
            
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
                            check = filt .* m;
                            check(check == 0) = [];

                            if size(check,2) < minimumPixels % Not enough pixels to continue
                                error('Not enough pixels in region. Please change mask radii.')
                            end

                            % Apply filter across full trajectory.
                            for n = 1:size(traj1,4)

                                if obj.meanDisksLogical(1,a) == 1 % Disk belongs

                                    S1 = traj1(:,:,1,n) .* filt; % reference img
                                    S1(S1 == 0) = [];

                                    q(a,b,n) = mean(S1(:)) / 255;
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