% Places local spatial discontinuities over wedge edges
classdef localDiscontinuities < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        % Basic parameters
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms

        % Regular gratings
        wedges = 2; % Number of wedges in baseline stimulus (2 = split-field spot)
        rotation = 0; % Rotate wedges as needed
        discontinuityDiameter = 50; % in microns. set to 0 to autosize to RF
        discontinuityLocation = [0 25 50 75]; % Location of spatial discontinuities (% of RF)

        % Stimulus parameters
        centerDiameter = 300; % in um
        contrast = 0.9; % 0 - 1
        temporalFrequency = 4 % Hz. Set to 0 for spatial flash
        backgroundIntensity = 0.168 % (0-1)

        % Additional
        randomizeOrder = true;
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        masks
        discontinuityMasks
        sequence
        counter
        centerDisk
        surroundDisk
        experimentTracker
        centerTracker
    end
    
    properties (Hidden, Transient)
        analysisFigure
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

            % Starting parameters
            centerRadius_px = obj.rig.getDevice('Stage').um2pix(obj.centerDiameter)/2;
            canvasSize      = obj.rig.getDevice('Stage').getCanvasSize();

            % Define polar space
            [xx,yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            r = sqrt((xx - canvasSize(1)/2).^2 + (yy - canvasSize(2)/2).^2);
            th = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2));
            th = abs(th-pi/2);
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);

            % Apply rotation
            th = mod(th + obj.rotation,360);
            
            % Define possible locations
            obj.centerDisk = (r <= centerRadius_px);

            % Build wedges
            obj.masks{1,1} = zeros(size(r));
            obj.masks{1,2} = zeros(size(r));
            A = obj.wedges;
            A(A == 0) = 1;
            theta = 0:360/A:360;
            for a = 1:length(theta) - 1
                m = (th >= theta(a) & th < theta(a+1));
                obj.masks{1,mod(a,2)+1} = obj.masks{1,mod(a,2)+1} + m;
            end

            % Build discontinuities
            obj.experimentTracker = [];
            for a = 1:length(theta)-1
                for b = 1:length(obj.discontinuityLocation)
                    if obj.discontinuityLocation(b) > 0 || (obj.discontinuityLocation(b) == 0 && a == 1)
                        rotatedTheta = mod(theta(a) + obj.rotation,360);
                        obj.experimentTracker = [obj.experimentTracker; rotatedTheta obj.discontinuityLocation(b)];
                    end
                end
            end
            
            % Add control
            obj.experimentTracker = [NaN NaN; obj.experimentTracker];

            estTime = size(obj.experimentTracker,1) * obj.numberOfAverages * ((obj.stimTime + obj.preTime + obj.tailTime + 500) / 1000 / 60);
            disp(strcat('estimated stimulus time:', mat2str(estTime),' min'))

            % Randomize order
            obj.sequence = 1:size(obj.experimentTracker,1);
            if obj.randomizeOrder == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end

            obj.counter = 0;
            obj.sequence = repmat(obj.sequence,1,obj.numberOfAverages);
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);

            epoch.addParameter('specificTheta', obj.experimentTracker(obj.sequence(obj.counter+1),1));
            epoch.addParameter('specificDistance', obj.experimentTracker(obj.sequence(obj.counter+1),2));
        end

        function p = createPresentation(obj)
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            centerRadius_px = obj.rig.getDevice('Stage').um2pix(obj.centerDiameter)/2;
            
            if obj.discontinuityDiameter > 0
                discontinuityRadius_px = obj.rig.getDevice('Stage').um2pix(obj.discontinuityDiameter)/2;
            else
                adjDiameter = obj.centerDiameter / length(obj.discontinuityLocation);
                discontinuityRadius_px = obj.rig.getDevice('Stage').um2pix(adjDiameter)/2;
            end
            
            disp(obj.experimentTracker(obj.sequence(obj.counter+1),:))
            theta = 180 - obj.experimentTracker(obj.sequence(obj.counter+1),1); % flip x-axis for monitor
            distance = obj.experimentTracker(obj.sequence(obj.counter+1),2);
            
            % Convert to euclidean space
            [xx,yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            rx = centerRadius_px .* distance/100 .* cos(deg2rad(theta));
            ry = centerRadius_px .* distance/100 .* sin(deg2rad(theta));
            r_tmp = sqrt((xx - canvasSize(1)/2 + rx).^2 + (yy - canvasSize(2)/2 + ry).^2);
            discontinuity = (r_tmp < discontinuityRadius_px);
            
            m = cat(3,obj.masks{1,1} .* (discontinuity == 0) + discontinuity .* obj.masks{1,2},...
                      obj.masks{1,2} .* (discontinuity == 0) + discontinuity .* obj.masks{1,1});

            % Add contrast reversing gratings
            for a = 1:2
                grate           = stage.builtin.stimuli.Grating('square'); % Square wave grating
                grate.size      = canvasSize;
                grate.position  = canvasSize / 2;
                grate.spatialFreq = 1e-5;
                grate.color     = 2 * obj.backgroundIntensity; % Amplitude of square wave
                grate.contrast  = obj.contrast; % Multiplier on square wave
                grateShape      = uint8(m(:,:,a) .* obj.centerDisk .* 255);
                grateMask       = stage.core.Mask(grateShape);
                grate.contrast  = obj.contrast;
                grate.setMask(grateMask);
              
                % Contrasting gratings between each set of masks
                if a == 1
                    grate.phase = 180;
                else
                    grate.phase = 0;
                end

                p.addStimulus(grate); % Add grating to the presentation
                
                % Control contrast
                if obj.temporalFrequency > 0
                    grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                        @(state)getGrateContrast(obj, state.time - obj.preTime/1e3));
                    p.addController(grateContrast); % Add the controller
                end
                
                % Hide during pre & post
                grateVisible = stage.builtin.controllers.PropertyController(grate, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(grateVisible);
            end
            obj.counter = obj.counter + 1;

            function c = getGrateContrast(obj, time)
                c = obj.contrast.*sin(2 * pi * obj.temporalFrequency * time);
            end
        end
 
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.sequence);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.sequence);
        end
        
    end
    
end