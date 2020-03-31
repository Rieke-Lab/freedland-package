% Uses weights from RFDiskArray to measure cell RF in relevant units
classdef RFWeighting < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        preTime = 250 % ms
        stimTime = 250 % ms
        tailTime = 250 % ms
        cellClass = 'ON' % type of cell
        centerSigma = [10,20,30,40,50,60,70,80,90,100] % in um
        annulusSize = [30,40,50,60,70,80,90,100] % in um
        maxLightLevel = 0.168; % maximum light intensity we encounter (0-1)
        randomizeOrder = true
        onlineAnalysis = 'none'
        numberOfAverages = uint16(2) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON', 'OFF'})        
        centerSigmaSequence
        annulusSizeSequence
        centerSigmaSelectionIndex
        annulusSizeSelectionIndex
        currentCenterSigma
        currentAnnulusSize
        backgroundIntensity
        spotIntensity
        imageMatrix
        trials
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
            
            if length(obj.centerSigma) > 1 && length(obj.annulusSize) > 1
                error('Please choose ranges for either "center sigma" or "annulus size", not both')
            end
            
            obj.centerSigmaSequence = obj.centerSigma;
            obj.annulusSizeSequence = obj.annulusSize;

            % Define x axis for online analysis
            if length(obj.centerSigmaSequence) > 1
                obj.showFigure('edu.washington.riekelab.freedland.figures.RFWeightingFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','center sigma');
                obj.trials = length(obj.centerSigmaSequence);
            else
                obj.showFigure('edu.washington.riekelab.freedland.figures.RFWeightingFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','annulus size');
                obj.trials = length(obj.annulusSizeSequence);
            end
            
            if obj.randomizeOrder == true
                obj.centerSigmaSequence = obj.centerSigmaSequence(randperm(length(obj.centerSigmaSequence)));
                obj.annulusSizeSequence = obj.annulusSizeSequence(randperm(length(obj.annulusSizeSequence)));
            end
            
            obj.centerSigmaSelectionIndex = 0;
            obj.annulusSizeSelectionIndex = 0;
            obj.currentCenterSigma = obj.centerSigmaSequence(obj.centerSigmaSelectionIndex+1);
            obj.currentAnnulusSize = obj.annulusSizeSequence(obj.annulusSizeSelectionIndex+1);
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);

            epoch.addParameter('specificCenterSigma', obj.currentCenterSigma);
            epoch.addParameter('specificAnnulusSize', obj.currentAnnulusSize);
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
            
            obj.centerSigmaSelectionIndex = mod(obj.centerSigmaSelectionIndex + 1,length(obj.centerSigma));
            obj.annulusSizeSelectionIndex = mod(obj.annulusSizeSelectionIndex + 1,length(obj.annulusSize));

            obj.currentCenterSigma = obj.centerSigmaSequence(obj.centerSigmaSelectionIndex+1);
            obj.currentAnnulusSize = obj.annulusSizeSequence(obj.annulusSizeSelectionIndex+1);
        end

        function lightLevelFilter = calculateFilter(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % Identify screen size

            % Convert to pixels
            centerSigmaPix = obj.currentCenterSigma ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            surroundSigmaPix = (obj.currentCenterSigma + obj.currentAnnulusSize) ./ obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');

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
                obj.backgroundIntensity = obj.maxLightLevel;
            else
                obj.backgroundIntensity = 0;
            end
        end
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages*obj.trials);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages*obj.trials);
        end
        
    end
    
end