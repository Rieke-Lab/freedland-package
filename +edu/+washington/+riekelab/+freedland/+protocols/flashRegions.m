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
        centerCuts   = 8;    % number of slices to divide center into
        rotate       = 0;    % degrees to rotate slices
        randomize    = true; % randomize order to present slices
        
        % Brightness
        diskIntensity = 0.319;       % 0 to 1
        backgroundIntensity = 0.168; % 0 to 1
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(3) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
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
            
            % Convert units
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            centerRadiusPix = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(...
                obj.centerRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'UM2PIX');
            
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
            rotations = 0:360/obj.centerCuts:360;
            obj.disks = zeros(size(r,1),size(r,2),obj.centerCuts);
            totalCombinations = 0;
            for a = 1:obj.centerCuts
                obj.disks(:,:,a) = (r <= centerRadiusPix) .* (th >= rotations(a) & th < rotations(a+1)) .* obj.diskIntensity;
                totalCombinations = totalCombinations + nchoosek(obj.centerCuts,a);
            end

            disp(strcat('total disk combinations: ',mat2str(totalCombinations)));
            totalTime = (totalCombinations * obj.numberOfAverages) * ((obj.preTime + obj.stimTime + obj.tailTime)*1.33/1000);
            disp(strcat('approx stimulus time (+33% rig delay):',mat2str(round(totalTime/60)),' minutes'));

            % Define all possible region combinations
            obj.selections = zeros(totalCombinations,obj.centerCuts);
            counter1 = 1;
            for a = 1:obj.centerCuts
                A = nchoosek(1:obj.centerCuts,a);
                for b = 1:size(A,1)
                    obj.selections(counter1,A(b,:)) = 1;
                    counter1 = counter1+1;
                end
            end

            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(totalCombinations);
            else
                obj.order = 1:totalCombinations;
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            % Flip for disks to correspond from 0 deg --> 360 deg
            % (counterclockwise from 3 o'clock)
            epoch.addParameter('flashedRegions',fliplr(obj.selections(obj.order(obj.counter+1),:)));
            
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
            specificDisks = logical(obj.selections(obj.order(obj.counter+1),:));
            image = sum(obj.disks(:,:,specificDisks),3);
            image(image == 0) = obj.backgroundIntensity;
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(uint8(image.*255));
            scene.size = [canvasSize(2) canvasSize(1)];
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