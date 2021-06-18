% Flashes various luminances in center and surround
classdef centerSurroundFlashes < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime     = 250 % ms
        stimTime    = 250 % ms
        tailTime    = 250 % ms
        
        % Cell information
        centerRadius     = 100; % radius of receptive-field center, in um
        centerContrast   = -90:10:90 % percent contrast
        surroundContrast = -90:10:90 % percent contrast
        backgroundIntensity = 0.168; % light intensity between trials (0-1)
        
        % Additional settings
        randomizeOrder = true
        onlineAnalysis = 'none'
        numberOfAverages = uint16(3) % number of repeats
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        centerSequence
        surroundSequence
        centerIndex
        surroundIndex
        specificCenter
        specificSurround
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

            % Sample all combinations of luminance values
            obj.centerSequence = repelem(obj.centerContrast,1,length(obj.surroundContrast));
            obj.surroundSequence = repmat(obj.surroundContrast,1,length(obj.centerContrast));
            if obj.randomizeOrder == true
                order = randperm(length(obj.centerSequence));
                obj.centerSequence = obj.centerSequence(order);
                obj.surroundSequence = obj.surroundSequence(order);
            end
            
            obj.centerIndex = 0;
            obj.surroundIndex = 0;
            obj.specificCenter = obj.centerSequence(obj.centerIndex + 1);
            obj.specificSurround = obj.surroundSequence(obj.surroundIndex + 1);
            
            % Estimate stimulus timing
            t = round(length(obj.centerSequence) .* (obj.preTime + obj.stimTime + obj.tailTime + 750)/1000 .* obj.numberOfAverages / 60,2); % +750ms rig delay
            disp(strcat('Est. stimulus time: ',mat2str(t), ' min.'))
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            % Save parameters as metadata
            epoch.addParameter('specificCenterContrast', obj.centerSequence(obj.centerIndex + 1));
            epoch.addParameter('specificSurroundContrast', obj.surroundSequence(obj.surroundIndex + 1));
        end

        function p = createPresentation(obj)

            % Stage parameters
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            p.setBackgroundColor(obj.backgroundIntensity);

            % Define settings
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            spotRadius = obj.rig.getDevice('Stage').um2pix(obj.centerRadius); % Convert from um to pixels
            
            obj.specificCenter = obj.centerSequence(obj.centerIndex + 1);
            obj.specificSurround = obj.surroundSequence(obj.surroundIndex + 1);
            centerLuminance = obj.backgroundIntensity .* (obj.specificCenter + 1);
            surroundLuminance = obj.backgroundIntensity .* (obj.specificSurround + 1);

            % Create surround            
            surround = stage.builtin.stimuli.Ellipse();
            surround.color = surroundLuminance;
            surround.radiusX = max(canvasSize);
            surround.radiusY = max(canvasSize);
            surround.position = canvasSize/2;
            p.addStimulus(surround);
            
            % Create center (on top of surround)            
            center = stage.builtin.stimuli.Ellipse();
            center.color = centerLuminance;
            center.radiusX = spotRadius;
            center.radiusY = spotRadius;
            center.position = canvasSize/2;
            p.addStimulus(center);
            
            % Add controllers to prevent flashing too early
            surroundVisible = stage.builtin.controllers.PropertyController(surround, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(surroundVisible);
            
            centerVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(centerVisible);
            
            % Adjust trackers
            obj.centerIndex = mod(obj.centerIndex + 1,length(obj.centerSequence));
            obj.surroundIndex = mod(obj.surroundIndex + 1,length(obj.surroundSequence));
        end    
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages * length(obj.centerSequence));
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages * length(obj.centerSequence));
        end
        
    end
    
end