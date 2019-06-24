% Provides a spinning oval stimulus that build a radial receptive field map
% simulatenously.
classdef RFRadialMap < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        stimTime = 3; % in sec. adjusts the growth rate.
        spinningFrequency = 1; % full rotations/sec. adjusts the spinning rate.
        ovalWidth = 30; % in um. fixed.
        ovalIntensity = 0.5; % light intensity of the spinning oval (0-1)
        backgroundIntensity = 0; % light intensity of background (0-1)
        onlineAnalysis = 'none'
        numberOfAverages = uint16(10) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        expandingSpeedPx
        timeTraj
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end
        
        function prepareRun(obj)

            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size
            canvasLimit = min(canvasSize) / 2; % stop at nearest edge of monitor.
            obj.expandingSpeedPx = canvasLimit / obj.stimTime; % in px/sec   
            
            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            obj.showFigure('edu.washington.riekelab.freedland.figures.RFRadialMapFigure',...
                obj.rig.getDevice(obj.amp),obj.rig.getDevice('Stage'),...
                obj.rig.getDevice('Frame Monitor'),'spinningFrequency',obj.spinningFrequency,...
                'expandingFrequency',obj.expandingSpeedPx,'canvasSize',canvasSize,...
                'ovalWidth',obj.rig.getDevice('Stage').um2pix(obj.ovalWidth),...
                'umPerPix',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));                
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = obj.stimTime;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
        end

        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size
            
            % identify stage parameters
            p = stage.core.Presentation(obj.stimTime);
            
            % calculate background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            % Build vectors for animation
            obj.timeTraj = 0:0.001:obj.stimTime; % 1 ms resolution
            
            spinningVec = mod(obj.timeTraj .* (obj.spinningFrequency * 360),360); % radial map
            expandingVec = (obj.timeTraj .* obj.expandingSpeedPx) + obj.rig.getDevice('Stage').um2pix(obj.ovalWidth); % RF map, starts at 2px
            
            % Place ellipse mask
            ellipse = stage.builtin.stimuli.Ellipse();
            ellipse.position = canvasSize/2; % center of monitor.
            ellipse.radiusY = obj.rig.getDevice('Stage').um2pix(obj.ovalWidth); % fixed width
            ellipse.color = obj.ovalIntensity;
            p.addStimulus(ellipse); % add mask
            
            sceneVisible = stage.builtin.controllers.PropertyController(ellipse, 'visible', ...
                @(state)state.time < (obj.stimTime*0.99)); % stop right before end to prevent adaptation
            p.addController(sceneVisible);
            
            % Animate spin
            ellipseSpin = stage.builtin.controllers.PropertyController(ellipse,...
                'orientation', @(state)spinSpinSpin(obj, state.time, spinningVec));
            p.addController(ellipseSpin); % spin mask
            
            % Animate growth
            ellipseGrow = stage.builtin.controllers.PropertyController(ellipse,...
                'radiusX', @(state)growGrowGrow(obj, state.time, expandingVec));
            p.addController(ellipseGrow); % grow mask
            
            function s = spinSpinSpin(obj, time, spinningVec)
                if time < 0 % do nothing
                    s = 0;
                elseif time > obj.timeTraj(end)
                    s = spinningVec(end); % hold last
                else
                    s = interp1(obj.timeTraj,spinningVec,time);
                end
            end
            
            function r = growGrowGrow(obj, time, expandingVec)
                if time < 0
                    r = obj.rig.getDevice('Stage').um2pix(obj.ovalWidth); % must have a starting length
                elseif time > obj.timeTraj(end)
                    r = expandingVec(end); % hold last
                else
                    r = interp1(obj.timeTraj,expandingVec,time);
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