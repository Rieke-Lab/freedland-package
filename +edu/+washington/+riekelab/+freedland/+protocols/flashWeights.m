% Use weights from flashRegions to design new spatial stimuli that elicit the SAME
% response as a uniform disk.
% By J. Freedland, 2020.
classdef flashWeights < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime  = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Cell type
        cellClass = 'ON'; % on cell or off cell?
        weights = [0.5 0.6 0.7 0.75 0.8 0.85 0.9 1]; % integration weight of each region (see: flashRegions)
        
        % Image information
        centerRadius = 100;  % in um
        cuts         = 8;    % number of slices to divide into
        cutLocation  = 'surround-only'; % locations to place cuts
        includeNegativeContrasts = true; % add negative contrast cases to stimulus
        rotate       = 0;    % degrees to rotate slices
        randomize    = true; % randomize order to present slices
        
        % Brightness
        contrast = 0.25; % 0 to 1
        backgroundIntensity = 0.168; % 0 to 1
        
        % Additional parameters
        numberOfStimuli = 10;
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(3) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON', 'OFF'}) 
        cutLocationType = symphonyui.core.PropertyType('char', 'row', {'center-only', 'surround-only','full-field'}) 
        disks
        selections
        order
        counter
        intensities
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
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis); 
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            obj.showFigure('edu.washington.riekelab.freedland.figures.receptiveFieldFitFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','experimentID');
            
            % Convert units
            canvasSize = fliplr(obj.rig.getDevice('Stage').getCanvasSize());
            centerRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.centerRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            
            %%% Create pixel space
            [xx,yy] = meshgrid(1:canvasSize(2),1:canvasSize(1));
            r = sqrt((xx - canvasSize(2)/2).^2 + (yy - canvasSize(1)/2).^2); 
            th = atan((xx - canvasSize(2)/2) ./ (yy - canvasSize(1)/2));
            th = abs(th-pi/2);              
    
            % Adjust theta space for strange monitors
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);
            th = mod(th + obj.rotate,360); % Rotate as required
            %%%

            % Build masks
            rotations = 0:360/obj.cuts:360;
            obj.disks = zeros(size(r,1),size(r,2),obj.cuts);
            if strcmp(obj.cutLocation,'center-only') || strcmp(obj.cutLocation,'full-field')
                m = (r <= centerRadiusPix);
            elseif strcmp(obj.cutLocation,'surround-only')
                m = (r > centerRadiusPix) & (r <= max(canvasSize));
            end
            
            % Split into regions
            for a = 1:obj.cuts
                obj.disks(:,:,a) = m .* (th >= rotations(a) & th < rotations(a+1));
            end
            
            % Adjust for special cases
            if strcmp(obj.cutLocation,'surround-only') % Add center disk to mix
                obj.disks = cat(3,(r <= centerRadiusPix),obj.disks);
            elseif strcmp(obj.cutLocation,'full-field') % Repeat slices for surround disks
                m = (r > centerRadiusPix) & (r <= max(canvasSize));
                extraDisks = zeros(size(r,1),size(r,2),obj.cuts);
                for a = 1:obj.cuts
                    extraDisks(:,:,a) = m .* (th >= rotations(a) & th < rotations(a+1));
                end
                obj.disks = cat(3,obj.disks,extraDisks);
            elseif strcmp(obj.cutLocation,'center-only') && length(obj.weights) == obj.cuts+1
                % Include a surround weight in the mix
                obj.disks = cat(3,obj.disks,(r > centerRadiusPix) & (r <= max(canvasSize)));
            end
            
            % Calculate intensities parameters
            obj.intensities = prepIntensities(obj);

            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(size(obj.intensities,1));
            else
                obj.order = 1:size(obj.intensities,1);
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('regionIntensity',(obj.intensities(obj.order(obj.counter+1),:)));
            epoch.addParameter('experimentID',obj.order(obj.counter+1));
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();     
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);

            % Identify specific parameters
            p.setBackgroundColor(obj.backgroundIntensity)  
            specificIntensity = obj.intensities(obj.order(obj.counter+1),:);
            
            % Build image
            image = zeros(size(obj.disks,1),size(obj.disks,2));
            for a = 1:obj.cuts
                image = image + obj.disks(:,:,a) .* specificIntensity(a);
            end
            image(sum(obj.disks,3) == 0) = obj.backgroundIntensity;
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(uint8(image.*255));
            scene.size = canvasSize;
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            % Only allow image to be visible during specific time
            p.addStimulus(scene);
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end
        
        function intensity = prepIntensities(obj)
            
            x = obj.weights';
            if strcmp(obj.cutLocation,'surround-only')
                x = x(2:end); % ignore center weight
            elseif strcmp(obj.cutLocation,'center-only') && length(obj.weights) == obj.cuts+1
                x = x(1:end-1); % ignore surround weight
            end
            B = sum(obj.contrast*x);
            
            % A * x = B
            % We have cell weights (x) and total spike counts (B). There are
            % infinite possibilities for A. So, we randomly generate a
            % number of givens and solve for the remaining values.
            remainingValues = 1;
            givens = length(x) - remainingValues;
            
            % Randomly generate contrasts
            % ~61% of cases will be unusable at 50% contrast
            % ~97% of cases will be unusable at 25%/75% contrast
            A = rand(round(obj.numberOfStimuli*1000),givens); % Enough for 0.1% usability
            
            % Solve for remaining contrasts
            A = [A B - A*x(1:givens) / x(givens+1:end)];
            A(sum(A<0,2)>0 | sum(A>1,2)>0,:) = []; % Ensure contrasts are between 0-1
            
            % Add static behavior to center/surround as needed
            if strcmp(obj.cutLocation,'surround-only')
                A = [repelem(obj.contrast,size(A,1),1), A];
            elseif strcmp(obj.cutLocation,'center-only') && length(obj.weights) == obj.cuts+1
                A = [A, repelem(obj.contrast,size(A,1),1)];
            end
            
            % Add case with uniform disk. If weights are correct, all
            % collections of intensities should elicit the same response.
            A = [repelem(obj.contrast,1,size(obj.masks,3));A(1:obj.numberOfStimuli-1,:)];
            
            % Sanity check
            % round(A * x,2) == round(B,2);
            
            if obj.includeNegativeContrasts == true
                intensity = cat(1,(A+1) .* obj.backgroundIntensity,(1-A) .* obj.backgroundIntensity);
            else
                if strcmp(obj.cellClass,'ON')
                    intensity = (A+1) .* obj.backgroundIntensity;
                elseif strcmp(obj.cellClass,'OFF')
                    intensity = (1-A) .* obj.backgroundIntensity;
                end
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.order);
        end
    end
end