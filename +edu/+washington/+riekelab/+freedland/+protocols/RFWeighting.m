% Uses weights from RFDiskArray to measure cell RF in relevant units
classdef RFWeighting < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        preTime = 250 % ms
        stimTime = 250 % ms
        tailTime = 250 % ms
        cellClass = 'ON' % type of cell
        centerSigma = [40 50 60 65 70 75 80 85 90 95 100 110 120] % in um
        annulusSize = 100 % in um
        maxLightLevel = 0.25; % maximum light intensity we encounter (0-1)
        randomizeOrder = false
        onlineAnalysis = 'none'
        numberOfAverages = uint16(2) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON', 'OFF'})        
        centerSigmaSequence
        selectionIndex
        currentCenterSigma
        backgroundIntensity
        spotIntensity
        imageMatrix
        centerRadius
        surroundRadius
        annulusRadius
    end
    
    properties (Hidden, Transient)
        
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
            obj.showFigure('edu.washington.riekelab.freedland.figures.RFWeightingFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,...
                'annulusSize',obj.annulusSize);

            % Create spot size sequence.
            obj.centerSigmaSequence = obj.centerSigma;
            
            if obj.randomizeOrder == true
                obj.centerSigmaSequence = obj.centerSigmaSequence(randperm(size(obj.centerSigmaSequence,2))); % place in random order
            end
            
            obj.selectionIndex = 1; % walks along random order
            obj.currentCenterSigma = obj.centerSigmaSequence(obj.selectionIndex);
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('centerSigma', obj.currentCenterSigma);
            epoch.addParameter('surroundSigma', obj.currentCenterSigma + obj.annulusSize);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity)
        end

        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % identify screen size
            
            % identify stage parameters
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % calculate background intensity
            rfFilter = calculateFilter(obj);
            p.setBackgroundColor(obj.backgroundIntensity);
            
            % scale to maximum light intensity
            obj.imageMatrix = uint8(rfFilter .* 255);
            
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
            
            obj.selectionIndex = obj.selectionIndex + 1; % walks along random order
            
            if obj.selectionIndex > size(obj.centerSigmaSequence,2)
                obj.selectionIndex = 1;
            end
            obj.currentCenterSigma = obj.centerSigmaSequence(obj.selectionIndex);
        end

        function lightLevelFilter = calculateFilter(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Identify screen size

            % Convert to pixels
            centerSigmaPix = obj.currentCenterSigma ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            surroundSigmaPix = (obj.currentCenterSigma + obj.annulusSize) ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');

            % Generate 2D gaussians
            centerGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],centerSigmaPix);
            surroundGaus = fspecial('gaussian',[canvasSize(2) canvasSize(1)],surroundSigmaPix);

            % Calculate difference of gaussians
            diffGaussian = centerGaus - surroundGaus;
            normFilter = diffGaussian ./ max(diffGaussian(:)); % Normalize filter

            % Shift gaussian s.t. all values are positive
            shift = abs(min(normFilter(:)));  
            normFilter = normFilter + shift;
            normFilter = normFilter ./ max(normFilter(:)); % Renormalize
            
            % This is the output from RFDiskArray.
            % However, we want to control light levels very closely.
            lightLevelFilter = normFilter .* obj.maxLightLevel;
            
            if strcmp(obj.cellClass,'OFF')
                lightLevelFilter = abs(lightLevelFilter - obj.maxLightLevel); % Complementary
            end
            
            obj.backgroundIntensity = mean(lightLevelFilter(:));
        end
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages*size(obj.centerSigmaSequence,2));
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages*size(obj.centerSigmaSequence,2));
        end
        
    end
    
end