classdef SingleSpotColorOpponency < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    % Presents a set of single spot stimuli to a Stage canvas and records from the specified amplifier.
    
    properties
        amp                             % Output amplifier
        preTime = 250                   % Spot leading duration (ms)
        stimTime = 250                  % Spot duration (ms)
        tailTime = 250                  % Spot trailing duration (ms)
        spotIntensity = 0.2             % Spot light intensity (0-1)
        spotDiameter = 300              % Spot diameter size (um)
        color1 = 'red';
        color2 = 'blue';
        backgroundIntensity = 0         % Background light intensity (0-1)
        numberOfAverages = uint16(5)    % Number of epochs
        interpulseInterval = 0          % Duration between spots (s)
    end
    
    properties (Hidden)
        ampType
        color1Type = symphonyui.core.PropertyType('char', 'row', {'red', 'green', 'blue'}) 
        color2Type = symphonyui.core.PropertyType('char', 'row', {'red', 'green', 'blue'}) 
    end
    
    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end
        
        function prepareRun(obj)
            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            
            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.figures.MeanResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.figures.FrameTimingFigure', obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            obj.showFigure('edu.washington.riekelab.freedland.figures.rfFlashFigureRatio',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,...
                'stimTime',obj.stimTime,'tailTime',obj.tailTime,'xval','figureLeader'); 
        end
        
        function p = createPresentation(obj)
            device = obj.rig.getDevice('Stage');
            canvasSize = device.getCanvasSize();
            
            spotDiameterPix = device.um2pix(obj.spotDiameter);
            
            p = stage.core.Presentation(2 * (obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            cycleTime = (obj.preTime + obj.stimTime + obj.tailTime) * 1e-3;
            p.setBackgroundColor(obj.backgroundIntensity);
            
            spot1 = stage.builtin.stimuli.Ellipse();
            spot1.radiusX = spotDiameterPix/2;
            spot1.radiusY = spotDiameterPix/2;
            spot1.position = canvasSize/2;
            
            if strcmp(obj.color1,'red')
                spot1.color = [obj.spotIntensity 0 0];
            elseif strcmp(obj.color1,'green')
                spot1.color = [0 obj.spotIntensity 0];
            elseif strcmp(obj.color1,'blue')
                spot1.color = [0 0 obj.spotIntensity];
            end

            p.addStimulus(spot1);
            
            spotVisible1 = stage.builtin.controllers.PropertyController(spot1, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(spotVisible1);
            
            spot2 = stage.builtin.stimuli.Ellipse();
            spot2.radiusX = spotDiameterPix/2;
            spot2.radiusY = spotDiameterPix/2;
            spot2.position = canvasSize/2;
            
            if strcmp(obj.color2,'red')
                spot2.color = [obj.spotIntensity 0 0];
            elseif strcmp(obj.color2,'green')
                spot2.color = [0 obj.spotIntensity 0];
            elseif strcmp(obj.color2,'blue')
                spot2.color = [0 0 obj.spotIntensity];
            end
            
            p.addStimulus(spot2);
            
            spotVisible2 = stage.builtin.controllers.PropertyController(spot2, 'visible', ...
                @(state)state.time >= cycleTime + obj.preTime * 1e-3 && state.time < cycleTime + (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(spotVisible2);
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = 2 * (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('figureLeader', 1);
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages;
        end
        
    end
    
end