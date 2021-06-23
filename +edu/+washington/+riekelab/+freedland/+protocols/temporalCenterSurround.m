% Flashes various luminances in center and surround with spatial contrast
% in center
classdef temporalCenterSurround < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime     = 250 % ms
        stepTime    = 250 % time between luminance steps (in ms)
        tailTime    = 250 % ms
        
        % Center-surround information
        centerRadius            = 100; % radius of receptive-field center, in um
        annulusSize             = 50;  % null region between center and surround, in um
        
        % Luminance information
        centerIntensity         = [-20 -10 0 10 20 30 40]; % luminance of uniform center, in percent contrast
        surroundIntensity       = [-90 -50 0 50 90]; % luminance of uniform surround, in percent contrast
        backgroundIntensity     = 0.168;             % light intensity between trials (0-1)
        
        % Additional settings
        randomizeOrder  = true
        onlineAnalysis  = 'extracellular'
        numberOfAverages = uint16(3) % number of repeats
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        centerLuminance
        surroundLuminance
        counter
        order
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
            
            % Define set of values
            obj.centerLuminance = repmat(obj.centerIntensity,1,length(obj.surroundIntensity));
            obj.surroundLuminance = repelem(obj.centerIntensity,1,length(obj.surroundIntensity));
            
            % Setup display
            obj.counter = 0;
            if obj.randomizeOrder == true
                obj.order = randperm(size(obj.centerLuminance,1));
            else
                obj.order = 1:size(obj.centerLuminance,1);
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);

            % Save parameters as metadata
            epoch.addParameter('specificCenter',obj.centerLuminance(obj.order(obj.counter+1)));
            epoch.addParameter('specificSurround',obj.surroundLuminance(obj.order(obj.counter+1)));
        end

        function p = createPresentation(obj)

            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();     
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity
            
            % Convert units
            spotRadiusPix = obj.rig.getDevice('Stage').um2pix(obj.centerRadius);
            annulusPix = obj.rig.getDevice('Stage').um2pix(obj.annulusSize);
            
            % Convert contrast to luminance
            spCenter = obj.backgroundIntensity .* (obj.centerLuminance(obj.order(obj.counter+1))/100 + 1);
            spSurround = obj.backgroundIntensity .* (obj.surroundLuminance(obj.order(obj.counter+1))/100 + 1);
            
            % Create center spot.            
            center = stage.builtin.stimuli.Ellipse();
            center.color = spCenter;
            center.radiusX = spotRadiusPix - annulusPix/2;
            center.radiusY = spotRadiusPix - annulusPix/2;
            center.position = canvasSize/2;
            p.addStimulus(center);
            
            % Create surround spot.
            surround = stage.builtin.stimuli.Rectangle();
            surround.position = canvasSize/2;
            surround.color = spSurround;
            surround.size = [max(canvasSize) max(canvasSize)];
            mask = stage.core.Mask.createAnnulus((aperatureDiameter-annulusPix/2)/max(canvasSize),...
                    (aperatureDiameter+annulusPix/2)/max(canvasSize), 1024);
            surround.setMask(mask);
            p.addStimulus(surround);

            %%% Temporal animation: 
                % 0 -> preTime: no stimulus
                % Step 1: Flash center
                % Step 2: Flash surround
                % Step 3: no stimulus
                % Step 4: Flash center & surround
                % Step 4 -> tailTime: no stimulus
            centerVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                (@(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stepTime.*2) * 1e-3) ||...
                (@(state)state.time >= (obj.preTime + obj.stepTime.*3) * 1e-3 && state.time < (obj.preTime + obj.stepTime.*4) * 1e-3));
            surroundVisible = stage.builtin.controllers.PropertyController(center, 'visible', ...
                (@(state)state.time >= (obj.preTime + obj.stepTime) * 1e-3 && state.time < (obj.preTime + obj.stepTime.*2) * 1e-3) ||...
                (@(state)state.time >= (obj.preTime + obj.stepTime.*3) * 1e-3 && state.time < (obj.preTime + obj.stepTime.*4) * 1e-3));
            
            p.addController(centerVisible);
            p.addController(surroundVisible);

            obj.counter = mod(obj.counter + 1,length(obj.order));
        end    
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages * length(obj.order));
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages * length(obj.order));
        end
        
    end
    
end