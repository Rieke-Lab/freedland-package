classdef contrastReversingGratingsCenterSurround < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms
        centerBarWidth = [5 10 20 40 60 80] % um
        centerContrast = [0.05 0.1 0.2 0.3 0.6] % relative to mean (0-1)
        surroundBarWidth = [40 80 120 160 200 240 300] % um
        surroundContrast = [0.05 0.1 0.2 0.3 0.6] % relative to mean (0-1)
        surroundPhase = 'in-phase';
        temporalFrequency = 4 % Hz
        apertureDiameter = 300; % um
        rotation = 0; % deg
        backgroundIntensity = 0.168 % (0-1)
        randomizeOrder = true;
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        surroundPhaseType = symphonyui.core.PropertyType('char', 'row', {'in-phase', 'out-of-phase'})
        sequence
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
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',1);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Arrange experiment conditions
            counter2 = 1;
            obj.sequence = zeros(length(obj.centerContrast).*length(obj.centerBarWidth).*length(obj.surroundContrast).*length(obj.surroundBarWidth),4);
            for a = 1:length(obj.centerContrast)
                for b = 1:length(obj.centerBarWidth)
                    for c = 1:length(obj.surroundContrast)
                        for d = 1:length(obj.surroundBarWidth)
                            obj.sequence(counter2,:) = [obj.centerContrast(a), obj.centerBarWidth(b), obj.surroundContrast(c), obj.surroundBarWidth(d)];
                            counter2 = counter2 + 1;
                        end
                    end
                end
            end
   
            % Randomize order
            if obj.randomizeOrder == true
                obj.sequence = obj.sequence(randperm(size(obj.sequence,1)),:);
            end

            obj.counter = 1;
            obj.sequence = repmat(obj.sequence,obj.numberOfAverages,1);
        end
             
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('currentCenterContrast', obj.sequence(obj.counter,1));
            epoch.addParameter('currentCenterBarWidth', obj.sequence(obj.counter,2));
            epoch.addParameter('currentSurroundContrast', obj.sequence(obj.counter,3));
            epoch.addParameter('currentSurroundBarWidth', obj.sequence(obj.counter,4));
        end

        function p = createPresentation(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            
            %convert from microns to pixels...
            apertureDiameterPix = obj.rig.getDevice('Stage').um2pix(obj.apertureDiameter);
            centerBarWidthPix   = obj.rig.getDevice('Stage').um2pix(obj.sequence(obj.counter,2));
            surroundBarWidthPix = obj.rig.getDevice('Stage').um2pix(obj.sequence(obj.counter,4));
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            
            % Create grating in RF center and surround
            for a = 1:2 % Center first, then surround
                grate  = stage.builtin.stimuli.Grating('square'); %square wave grating
                grate.orientation = obj.rotation;
                grate.position    = canvasSize/2;
                grate.color       = 2*obj.backgroundIntensity;
                
                if a == 1 % Center
                    grate.size = [apertureDiameterPix, apertureDiameterPix];
                    grate.spatialFreq = 1/(2*centerBarWidthPix); %convert from bar width to spatial freq
                elseif a == 2 % Surround
                    grate.size = [max(canvasSize) max(canvasSize)];
                    if surroundBarWidthPix > 0
                        grate.spatialFreq = 1/(2*surroundBarWidthPix); %convert from bar width to spatial freq
                    end
                    % Apply circular hole for center
                    mask = stage.core.Mask.createCircularAperture(apertureDiameterPix/max(canvasSize), 1024); %circular aperture
                    grate.setMask(mask);
                end

                %calc to apply phase shift s.t. a contrast-reversing boundary
                %is in the center regardless of spatial frequency. Arbitrarily
                %say boundary should be positve to right and negative to left
                %crosses x axis from neg to pos every period from 0
                zeroCrossings   = 0:(grate.spatialFreq^-1):grate.size(1); 
                offsets         = zeroCrossings-grate.size(1)/2; %difference between each zero crossing and center of texture, pixels
                [shiftPix, ~]   = min(offsets(offsets>0)); %positive shift in pixels
                phaseShift_rad  = (shiftPix/(grate.spatialFreq^-1))*(2*pi); %phaseshift in radians
                phaseShift      = 360*(phaseShift_rad)/(2*pi); %phaseshift in degrees
                grate.phase = phaseShift; 
                p.addStimulus(grate);

                % Make it contrast-reversing
                if (obj.temporalFrequency > 0) 
                    if a == 1 % Center
                        grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                            @(state)getGrateContrastCenter(obj, state.time - obj.preTime/1e3));
                    else
                        grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                            @(state)getGrateContrastSurround(obj, state.time - obj.preTime/1e3));
                    end
                    p.addController(grateContrast); %add the controller
                end
                
                grateVisible = stage.builtin.controllers.PropertyController(grate, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(grateVisible);
            end
            
            function c = getGrateContrastCenter(obj, time)
                c = obj.sequence(obj.counter,1).*sin(2 * pi * obj.temporalFrequency * time);
            end
            
            function c = getGrateContrastSurround(obj, time)
                if strcmp(obj.surroundPhase,'in-phase')
                    c = obj.sequence(obj.counter,3).*sin(2 * pi * obj.temporalFrequency * time);
                elseif strcmp(obj.surroundPhase,'out-of-phase')
                    c = obj.sequence(obj.counter,3).*sin(2 * pi * obj.temporalFrequency * time) .* -1;
                end
            end
            
            obj.counter = obj.counter + 1;
        end
 
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (size(obj.sequence,1)-1); 
            %%% To fix: last epoch doesn't start; quick fix is to simply
            %%% not play very last protocol (so protocol finishes)
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (size(obj.sequence,1)-1);
        end
        
        
    end
    
end