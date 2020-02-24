% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFWhiteNoiseDisk < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 500 % in ms
        tailTime = 250 % in ms
        
        % White noise information
        meanLightIntensity = 0.168; % average light intensity of noise
        maxContrast = 0.5; % maximum contrast for noise
        pixelSteps = 3; % how many discrete pixel values to include     
        
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
        disksIgnoreCut = [0 0]; % starting from the center disk and moving outwards, which disks should we NOT cut (keep circular)?

        % Disk type
        meanDisks = [1 2 3]; % starting from the center disk and moving outwards, which disks should be averaged?
        backgroundDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be left at background intensity?
        meanIntegration = 'gaussian'; % type of linear integration

        randomize = true; % plays sequence of protocols in random order
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        trajectoryType = symphonyui.core.PropertyType('char', 'row', {'natural','disk','both'})
        meanIntegrationType = symphonyui.core.PropertyType('char', 'row', {'uniform','gaussian'})
        overrideCoordinateType = symphonyui.core.PropertyType('char', 'row', {'pixels','RF'})
        backgroundIntensity
        imageDatabase
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
        frameTraj
        sequence
        selection
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)

            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            
            % Must have prerender set to avoid lagging through the replacement trajectory.     
%             if obj.rig.getDevice('Stage').getConfigurationSetting('prerender') == 0
%                 error('Must have prerender set') 
%             end
            
            % Must have at least 2 pixel types to allow consistent average
            if obj.pixelSteps < 2
                error('Please check pixelSteps. Insufficient number')
            end

            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));

            % Catch common errors.
            checkDisks = sum(1:obj.disks);
            checkAssignments = sum([obj.meanDisks obj.backgroundDisks]);
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
            
            % Estimate number of frames to generate
            monitorRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
            frames = round((obj.stimTime + 50) / 1000 * monitorRate);
                        
            % Build white noise movie
            img = whiteNoiseMovie(obj,1);
            img2 = img.*255;
            obj.imageDatabase = uint8(img2);  
            
            % Associate time with frame number
            obj.timeTraj = 0:(frames * monitorRate) - 1; 
            obj.frameTraj = repelem(1:frames,monitorRate);

            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Calculate screen size
            
            if strcmp(obj.overrideCoordinate,'RF') && obj.overrideRadiiLogical == true
                RFConversion(obj) % change coordinate systems
                obj.radii = obj.overrideRadii;
            else
                obj.radii = [0 cumsum(repelem(max(canvasSize)/(obj.disks),1,(obj.disks)))] ./ 2; % Evenly distributed (in VH)
            end
            
            % Recalculate masks in monitor pixels.
            obj.diskExpandRadius = 0.5;   % expands disks radially by percentage to ensure overlap.  
            obj.diskExpandAngle = 0; % expands disks' angle coverage by percentage to ensure overlap.  
            obj.diskOpacity = 1;    % opacity of disks placed. This is a good setting to play with during stimulus testing. (Matches pixels underneath).
                
            % Define space
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
            
            if obj.xSliceFrequency > 0 && obj.ySliceFrequency == 0
                thetaX = 90 / obj.xSliceFrequency;
                t = 0:thetaX:90-(1E-10);
                obj.theta = [t t + 180 360];
            elseif obj.xSliceFrequency == 0 && obj.ySliceFrequency > 0
                thetaY = 90 / obj.ySliceFrequency;
                t = 90:thetaY:180-(1E-10);
                obj.theta = [t t + 180 450];
                obj.theta(round(imgSize(1)/2):end,round(imgSize(2)/2):end) = obj.theta(round(imgSize(1)/2):end,round(imgSize(2)/2):end) + 360; % add
            elseif obj.xSliceFrequency > 0 && obj.ySliceFrequency > 0
                thetaX = 90 / obj.xSliceFrequency;
                thetaY = 90 / obj.ySliceFrequency;
                t = [0:thetaX:90-(1E-10) 90:thetaY:180-(1E-10)];
                obj.theta = [t t + 180 360];
            end

            % The case with no slices
            if sum([obj.xSliceFrequency obj.ySliceFrequency]) == 0

                obj.masks = zeros(canvasSize(2),canvasSize(1),size(obj.radii,2));
                obj.specificOpacity = 1;

                for a = 1:size(obj.radii,2) - 1
                    tempMask = m >= (obj.radii(1,a).*(1-obj.diskExpandRadius/100))...
                        & m <= (obj.radii(1,a+1).*(1+obj.diskExpandRadius/100));
                    obj.masks(:,:,a) = 1 - tempMask; 
                end

            else % The case with slices

                obj.masks = zeros(canvasSize(2),canvasSize(1),(size(obj.radii,2)-1).*(size(obj.theta,2) - 1));
                obj.specificOpacity = 1;
                counter = 1;
                
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
                        end

                        tempMask = dist .* th;
                        obj.masks(:,:,counter) = 1 - tempMask; 
                        counter = counter + 1;
                    end
                end              
            end
            
            obj.selection = 0;
            if obj.randomize == true
                obj.sequence = randperm(size(obj.masks,3) * size(obj.masks,4));
            else
                obj.sequence = 1: (size(obj.masks,3));
            end
        end
        
        function prepareEpoch(obj, epoch)
            
             prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('radii', obj.radii); % in pixels
            epoch.addParameter('diskNumber',obj.sequence(obj.selection + 1));
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);

            % Set background intensity
            p.setBackgroundColor(obj.meanLightIntensity);

            % Prep to display image
            scene = stage.builtin.stimuli.Image(obj.imageDatabase(:,:,1));
            scene.size = [size(obj.imageDatabase,2),size(obj.imageDatabase,1)];
            scene.position = canvasSize/2;

            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            
%             % Change static image appropriately
%             sceneImage = stage.builtin.controllers.PropertyController(scene,...
%                 'imageMatrix', @(state)getImage(obj, state.time - obj.preTime/1e3));
%             
%             function p = getImage(obj, time)
%                 if time <= 0
%                     p = obj.imageMatrix(:,:,1);
%                 elseif time > obj.timeTraj(end) % Hang on last frame
%                     p = obj.imageMatrix(:,:,end);
%                 else % Present new images as needed
%                     t = round(interp1(obj.timeTraj,obj.frameTraj,time));
%                     p = obj.imageMatrix(:,:,t);
%                 end
%             end
            
            p.addStimulus(scene);
%             p.addController(sceneImage);

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);

            %%%%%% Apply masks %%%%%% 
    
            annulus = stage.builtin.stimuli.Grating();
            annulus.position = canvasSize/2;
            annulus.color = obj.meanLightIntensity.*2;
            annulus.size = [canvasSize(1) canvasSize(2)]; 
            annulus.orientation = 0; 
            annulus.spatialFreq = 0;
            annulus.opacity = obj.specificOpacity;
            annulusA = uint8(obj.masks(:,:,obj.sequence(obj.selection+1)).*255);
            annulusMaskX = stage.core.Mask(annulusA);
            annulus.setMask(annulusMaskX);
            
            p.addStimulus(annulus);
            sceneVisible = stage.builtin.controllers.PropertyController(annulus, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            obj.selection = mod(obj.selection + 1,size(obj.masks,3));
        end

        function mov = whiteNoiseMovie(obj,frames)
            
            % Spectrum of possible pixel values
            pixelRange = (obj.meanLightIntensity - (obj.maxContrast * obj.meanLightIntensity)) : 1e-3 : (obj.meanLightIntensity + (obj.maxContrast * obj.meanLightIntensity));
            
            % Pull corresponding number of steps
            if obj.pixelSteps == 2
                pixelVals = pixelRange( [1, length(pixelRange)] );
            else
                pixelVals = pixelRange( [1, round(length(pixelRange) / (obj.pixelSteps - 1)) : round(length(pixelRange) / (obj.pixelSteps - 1)) : length(pixelRange)-1, length(pixelRange)] );
            end
            
            randVals = 0 : 1/(obj.pixelSteps) : 1; % Bin values
            
            % Make image
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            mov = zeros(canvasSize(2),canvasSize(1),frames);
                
            for b = 1:frames
                startingImg = rand(canvasSize(2),canvasSize(1));
                for a = 1:length(pixelVals)
                    mov(:,:,b) = mov(:,:,b) + double(startingImg <= randVals(a+1) & startingImg >= randVals(a)) .* pixelVals(a);
                end
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
            
            centerAnnulusBorder = min(sizing);
            annulusSurroundBorder = max(sizing);
            
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
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages .* size(obj.masks,3);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages .* size(obj.masks,3);
        end
    end
end