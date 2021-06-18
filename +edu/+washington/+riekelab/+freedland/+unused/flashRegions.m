% Flash regions in low-dimensional space to find correlations.
% By J. Freedland, 2019.
classdef flashRegions < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime  = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % Image information
        centerRadius = 100;  % in um
        slices       = 8;    % number of spatial slices to divide into
        sliceLocation = 'center'; % location of flashing slices: center or surround
        rotate       = 0;    % degrees to rotate slices
        randomize    = true; % randomize flashing order
        
        % Brightness
        diskIntensity       = 0.319; % 0 to 1
        backgroundIntensity = 0.168; % 0 to 1
        
        % Add-ons
        includeNegativeIntensity = true; % Add negative-contrast regions.
        includeSurroundIntensity = true; % For center-only stimuli, add surround flashes.
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(1) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        sliceLocationType = symphonyui.core.PropertyType('char', 'row', {'center', 'surround'}) 
        disks
        selections
        order
        counter
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
            obj.showFigure('edu.washington.riekelab.freedland.figures.flashRegionsFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','regions');
            
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
            rotations = 0:360/obj.slices:360;
            obj.disks = zeros(size(r,1),size(r,2),obj.slices);
            totalCombinations = 0;
            for a = 1:obj.slices
                if strcmp(obj.sliceLocation,'center')
                    obj.disks(:,:,a) = (r <= centerRadiusPix) .* (th >= rotations(a) & th < rotations(a+1));
                elseif strcmp(obj.sliceLocation,'surround')
                    obj.disks(:,:,a) = (r > centerRadiusPix & r <= min(canvasSize)/2) .* (th >= rotations(a) & th < rotations(a+1));
                end
                totalCombinations = totalCombinations + nchoosek(obj.slices,a);
            end
            
            % Add uniform region
            if strcmp(obj.sliceLocation,'center')
                obj.disks = cat(3,obj.disks,(r > centerRadiusPix & r <= min(canvasSize)/2));
            elseif strcmp(obj.sliceLocation,'surround')
                obj.disks = cat(3,(r <= centerRadiusPix),obj.disks);
            end

            % Define all possible combinations
            obj.selections = zeros(totalCombinations,obj.slices);
            counter1 = 1;
            for a = 1:obj.slices
                A = nchoosek(1:obj.slices,a);
                for b = 1:size(A,1)
                    obj.selections(counter1,A(b,:)) = 1;
                    counter1 = counter1+1;
                end
            end
            
            %%%%%%%%% For adding negative contrast central disks.
            % Adding a negative option for every combination of disks
            % produces an extremely large set. To minimize our experiment:
            negativeDiskDimension = 2; % Dimensional space to add negative disks for.
            if obj.includeNegativeIntensity == true
                % Identify all trials at proper dimension
                dimSelect = obj.selections(sum(obj.selections,2) == negativeDiskDimension,:);
                for a = 1:obj.slices
                    v = ones(size(dimSelect,1),size(dimSelect,2));
                    v(:,a) = -1; % Negate one disk
                    obj.selections = [obj.selections; dimSelect .* v];
                end
                obj.selections = unique(obj.selections,'rows');
            end

            %%%%%%%%% For adding disks outside of our specific region.
            if strcmp(obj.sliceLocation,'center')
                % Exclude surround for most of dataset
                obj.selections = [obj.selections,repelem(0,size(obj.selections,1),1)];
                if obj.includeSurroundIntensity == true
                    surroundDiskDimension = 2; % Add surround for specific dimension only
                    dimSelect = obj.selections(sum(obj.selections == -1,2) == 0 &...
                        sum(obj.selections,2) == surroundDiskDimension,:); 
                    dimSelect(:,end) = 1;
                    obj.selections = [obj.selections; dimSelect];
                end
            elseif strcmp(obj.sliceLocation,'surround')
                % Include center for entire experiment (helps measure inhibition)
                obj.selections = [repelem(1,size(obj.selections,1),1), obj.selections];
                
                % Add center-only case & surround-only case to keep proper rank
                centerOnly = [1 zeros(1,size(obj.selections,2)-1)];
                surroundOnly = [0 ones(1,size(obj.selections,2)-1)];
                obj.selections = [obj.selections; centerOnly; surroundOnly];
            end
            disp(strcat('total disk combinations: ',mat2str(size(obj.selections,1))));
            totalTime = (size(obj.selections,1) * obj.numberOfAverages) * ((obj.preTime + obj.stimTime + obj.tailTime + 1500)/1000);
            disp(strcat('approx stimulus time (+1.5 sec delay):',mat2str(round(totalTime/60)),' minutes'));

            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(size(obj.selections,1));
            else
                obj.order = 1:size(obj.selections,1);
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            epoch.addParameter('flashedRegions',obj.selections(obj.order(obj.counter+1),:));
            
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

            % Set image
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity
            
            % Image with positive contrast
            positiveDisks = obj.selections(obj.order(obj.counter+1),:) == 1;
            posImage = sum(obj.disks(:,:,positiveDisks),3) .* obj.diskIntensity;
            
            % Image with negative contrast
            contrast = abs(obj.backgroundIntensity - obj.diskIntensity);
            if obj.diskIntensity > obj.backgroundIntensity % ON Cell
                oppDiskIntensity = (obj.backgroundIntensity - contrast);
            else
                oppDiskIntensity = (obj.backgroundIntensity + contrast);
            end
            negativeDisks = obj.selections(obj.order(obj.counter+1),:) == -1;
            negImage = sum(obj.disks(:,:,negativeDisks),3) .* oppDiskIntensity;
            
            % Combine images
            image = posImage + negImage;
            image(image == 0) = obj.backgroundIntensity;
            
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
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.order);
        end
    end
end