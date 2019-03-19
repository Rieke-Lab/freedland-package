classdef RFEllipse < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        preTime = 250 % ms
        stimTime = 250 % ms
        tailTime = 250 % ms
        spotDistances = [40 60 80]; % distance from center (um)
        resolution = 15; % in degrees
        spotSize = 50; % um
        spotIntensity = 0.8; % maximum light intensity we encounter (0-1)
        backgroundIntensity = 0;
        compressVertically = 1; % form ellipse by compressing up/down direction. makes horizontal ellipse.
        compressHorizontally = 1; % form ellipse by compressing left/right direction. makes vertical ellipse.
        rotation = 0; % rotation of ellipse (counterclockwise) in deg.
        randomizeRotation = true % rotates randomly
        onlineAnalysis = 'none'
        numberOfAverages = uint16(1) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON', 'OFF'})   
        locationSeriesType = symphonyui.core.PropertyType('char', 'row', {'N S E W', 'NE NW SE SW', 'all'})  
        spotDistanceSequence
        selectionIndexD
        selectionIndexL
        spotLocation
        spotLocationSequence
        currentDistance
        currentLocation
        imageMatrix
        imageAnalysis
        holdover
        circSize
        circRadius
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
            obj.showFigure('edu.washington.riekelab.freedland.figures.RFEllipticalFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,...
                'spotDistances',obj.spotDistances,'resolution',obj.resolution);
            
            obj.spotDistanceSequence = obj.spotDistances;
            
            % we will assume a symmetric RF field (response at 0deg = 180deg).
            % this cuts our sampling in half!
            obj.spotLocationSequence = 0:obj.resolution:180-(1E-9); % in polar coordinates
                       
            if obj.randomizeRotation == true
                obj.spotLocationSequence = obj.spotLocationSequence(randperm(size(obj.spotLocationSequence,2))); % place in random order
            end
            
            obj.selectionIndexD = 1; % walks along random order
            obj.selectionIndexL = 1; 
            obj.currentDistance = obj.spotDistanceSequence(obj.selectionIndexD);
            obj.currentLocation = obj.spotLocationSequence(obj.selectionIndexL);
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            obj.imageAnalysis = zeros(canvasSize(2),canvasSize(1));
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('distance', obj.currentDistance);
            epoch.addParameter('angle', obj.currentLocation)
            obj.imageAnalysis = obj.imageAnalysis; % for figure
        end

        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size
            
            % identify stage parameters
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % calculate background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            % create circle for stimulus
            r = obj.rig.getDevice('Stage').um2pix(obj.spotSize);
            obj.circSize = r+1;
            obj.circRadius = r/2;
            circ = createCircle(obj);
            
            % find location
            theta = -deg2rad(obj.currentLocation); % switch for rotation
            d = obj.rig.getDevice('Stage').um2pix(obj.currentDistance);
            xx = round(sin(theta) * (d / obj.compressHorizontally) + (canvasSize(2) / 2));
            yy = round(cos(theta) * (d / obj.compressVertically) + (canvasSize(1) / 2));
            xx2 = (canvasSize(2) / 2) - round(sin(theta) * (d / obj.compressHorizontally));
            yy2 = (canvasSize(1) / 2) - round(cos(theta) * (d / obj.compressVertically));
            imgMatrix = ones(canvasSize(2),canvasSize(1))*obj.backgroundIntensity; % background intensity for img.
            imgMatrix(xx - obj.circRadius : xx + obj.circRadius,yy - obj.circRadius : yy + obj.circRadius) = circ;
            imgMatrix(xx2 - obj.circRadius : xx2 + obj.circRadius,yy2 - obj.circRadius : yy2 + obj.circRadius) = circ;
            
            obj.imageMatrix = uint8(imgMatrix*255); % scale to max light intensity of monitor.
            
            scene = stage.builtin.stimuli.Image(obj.imageMatrix);
            scene.size = [size(obj.imageMatrix,2) size(obj.imageMatrix,1)];
            scene.position = canvasSize/2;
            p.addStimulus(scene);
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            
            obj.selectionIndexL = obj.selectionIndexL + 1; % iterate by location
            
            if obj.selectionIndexL > size(obj.spotLocationSequence,2)
                obj.selectionIndexL = 1;
                obj.selectionIndexD = obj.selectionIndexD + 1; % then by radius
            end
            
            if obj.selectionIndexD > size(obj.spotDistanceSequence,2) % restart for next run
                obj.selectionIndexL = 1;
                obj.selectionIndexD = 1;
                if obj.randomizeRotation == true
                    obj.spotLocationSequence = obj.spotLocationSequence(randperm(size(obj.spotLocationSequence,2))); % place in random order
                end
            end
            
            obj.currentDistance = obj.spotDistanceSequence(obj.selectionIndexD);
            obj.currentLocation = obj.spotLocationSequence(obj.selectionIndexL);
        end
        
        function circ = createCircle(obj)

            [xx, yy] = meshgrid(1:obj.circSize,1:obj.circSize);
            m = sqrt((xx-(obj.circSize/2)).^2+(yy-(obj.circSize/2)).^2);
            
            circ = m <= obj.circRadius; % circ
            surr = m > obj.circRadius;
            circ = double(circ) .* obj.spotIntensity;
            surr = double(surr) .* obj.backgroundIntensity; % background in circ.
            circ = circ + surr;
        end
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages*size(obj.spotLocationSequence,2)*size(obj.spotDistanceSequence,2));
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages*size(obj.spotLocationSequence,2)*size(obj.spotDistanceSequence,2));
        end

    end
    
end