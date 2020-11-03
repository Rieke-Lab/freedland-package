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
        addNegativeIntensity = true; % Add trials with opposite-contrast adjacent disks
        
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
            obj.showFigure('edu.washington.riekelab.freedland.figures.flashRegionsFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','regions');
            
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
                obj.disks(:,:,a) = (r <= centerRadiusPix) .* (th >= rotations(a) & th < rotations(a+1));
                totalCombinations = totalCombinations + nchoosek(obj.centerCuts,a);
            end

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
            
            % Add opposing disk conditions
            if obj.addNegativeIntensity == true
                counter1 = 1;
                negativeIntensity = zeros(totalCombinations,obj.centerCuts);
                for a = 2 % DIMENSIONAL SPACE TO INCLUDE NEGATIVE DISKS IN
                    A = nchoosek(1:obj.centerCuts,a);
                    for b = 1:size(A,1) % Specific combination of disks
                        B = A(b,:);
                        for c = 1:a-1 % Find every possible combination of negative disks
                            C = nchoosek(1:a,c);
                            for d = 1:size(C,1)
                                negativeIntensity(counter1,B) = 1;
                                negativeIntensity(counter1,B(C(d,:))) = -1;
                                counter1 = counter1+1;
                            end
                        end
                    end
                end
                negativeIntensity(sum(abs(negativeIntensity),2) == 0,:) = [];
                obj.selections = [obj.selections; negativeIntensity];
                totalCombinations = size(obj.selections,1);
            end
            
            disp(strcat('total disk combinations: ',mat2str(totalCombinations)));
            totalTime = (totalCombinations * obj.numberOfAverages) * ((obj.preTime + obj.stimTime + obj.tailTime + 1500)/1000);
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
            positiveDisks = obj.selections(obj.order(obj.counter+1),:) == 1; % Already at intensity
            posImage = sum(obj.disks(:,:,positiveDisks),3) .* obj.diskIntensity;
            
            negativeDisks = obj.selections(obj.order(obj.counter+1),:) == -1;
            
            contrast = abs(obj.backgroundIntensity - obj.diskIntensity);
            if obj.diskIntensity > obj.backgroundIntensity % ON Cell
                oppDiskIntensity = (obj.backgroundIntensity - contrast);
            else
                oppDiskIntensity = (obj.backgroundIntensity + contrast);
            end
            negImage = sum(obj.disks(:,:,negativeDisks),3) .* oppDiskIntensity;
            
            image = posImage + negImage;
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