% Uses weights from RFDiskArray to measure cell RF in relevant units
classdef receptiveFieldFit < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        preTime = 250 % ms
        stimTime = 250 % ms
        tailTime = 250 % ms
        cellClass = 'ON' % type of cell
        centerSigma = [10,20,30,40,50,60,70,80,90,100,110,120] % in um
        annulusSize = [30,40,50,60,70,80,90,100,110,120,130,140,150] % in um
        backgroundIntensity = 0.168; % maximum light intensity we encounter (0-1)
        centerContrast = 0.5; % 0-1 for spot brightness
        surroundContrast = 0.5; % 0-1 for spot brightness
        disk = 'uniform';
        randomizeOrder = false
        onlineAnalysis = 'none'
        numberOfAverages = uint16(2) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON', 'OFF'})   
        diskType = symphonyui.core.PropertyType('char', 'row', {'uniform', 'spatial'})   
        centerSigmaSequence
        annulusSizeSequence
        centerSigmaSelectionIndex
        annulusSizeSelectionIndex
        rfSigmaCenter
        rfSigmaSurround
        spotIntensity
        imageMatrix
        trials
        micronsPerPixel
        monitorSize
        monitorFrameRate
        videoSize
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
                obj.showFigure('edu.washington.riekelab.freedland.figures.receptiveFieldFitFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','center sigma (um)');
                obj.trials = length(obj.centerSigmaSequence);
            else
                obj.showFigure('edu.washington.riekelab.freedland.figures.receptiveFieldFitFigure',...
                obj.rig.getDevice(obj.amp),'preTime',obj.preTime,'stimTime',obj.stimTime,'type','surround sigma (um)');
                obj.trials = length(obj.annulusSizeSequence);
            end
            
            if obj.randomizeOrder == true
                obj.centerSigmaSequence = obj.centerSigmaSequence(randperm(length(obj.centerSigmaSequence)));
                obj.annulusSizeSequence = obj.annulusSizeSequence(randperm(length(obj.annulusSizeSequence)));
            end
            
            obj.centerSigmaSelectionIndex = 1;
            obj.annulusSizeSelectionIndex = 1;
            obj.rfSigmaCenter = obj.centerSigmaSequence(obj.centerSigmaSelectionIndex);
            obj.rfSigmaSurround = obj.rfSigmaCenter + obj.annulusSizeSequence(obj.annulusSizeSelectionIndex);
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('rfSigmaCenter', obj.rfSigmaCenter);
            epoch.addParameter('rfSigmaSurround', obj.rfSigmaSurround);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity)
        end

        function p = createPresentation(obj)

            % identify stage parameters
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            p.setBackgroundColor(obj.backgroundIntensity);
            
            % Pull parameters for package
            obj.micronsPerPixel = obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.monitorSize = obj.rig.getDevice('Stage').getCanvasSize();
            obj.monitorSize = fliplr(obj.monitorSize); % Adjust to [height, width]
            obj.monitorFrameRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
            obj.videoSize = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(obj.monitorSize,obj.micronsPerPixel,'PIX2VH');

            % Pull filter
            rfFilter = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.calculateFilter(obj,false);

            % Scale excitatory/inhibitory regions accordingly.
            if strcmp(obj.cellClass,'ON')
                center = double(rfFilter > 0) .* (obj.backgroundIntensity .* (1+obj.centerContrast));
                surround = double(rfFilter <= 0) .* (obj.backgroundIntensity .* (1-obj.surroundContrast));
                
                % center = sum(center * RF) / sum(RF)
                if strcmp(obj.disk,'spatial')
                    center = obj.backgroundIntensity .* (1+obj.centerContrast) ./ rfFilter;
                    surround = obj.backgroundIntensity .* (1-obj.surroundContrast) ./ rfFilter;
                    center(center > 1 | center < 0) = 0;
                    surround(surround > 0 | surround < -1) = 0;
                    surround = abs(surround);
                end
                
                if obj.surroundContrast == 0
                    surround = double(rfFilter <= 0) .* (obj.backgroundIntensity .* (1-obj.surroundContrast));
                end
                
                if obj.centerContrast == 0
                    center = double(rfFilter > 0) .* (obj.backgroundIntensity .* (1+obj.centerContrast));
                end
                
            elseif strcmp(obj.cellClass,'OFF')
                center = double(rfFilter > 0) .* (obj.backgroundIntensity .* (1-obj.centerContrast));
                surround = double(rfFilter <= 0) .* (obj.backgroundIntensity .* (1+obj.surroundContrast));
                
                % center = sum(center * RF) / sum(RF)
                if strcmp(obj.disk,'spatial')
                    center = obj.backgroundIntensity .* (1-obj.centerContrast) ./ rfFilter;
                    surround = obj.backgroundIntensity .* (1+obj.surroundContrast) ./ rfFilter;
                    center(center > 1 | center < 0) = 0;
                    surround(surround > 0 | surround < -1) = 0;
                    surround = abs(surround);
                end
                
                if obj.surroundContrast == 0
                    surround = double(rfFilter <= 0) .* (obj.backgroundIntensity .* (1+obj.surroundContrast));
                end
                
                if obj.centerContrast == 0
                    center = double(rfFilter > 0) .* (obj.backgroundIntensity .* (1-obj.centerContrast));
                end
            end
            
            % Scale to maximum light intensity
            obj.imageMatrix = uint8((center+surround) .* 255);
            
            % Display in symphony
            scene = stage.builtin.stimuli.Image(obj.imageMatrix);
            scene.size = fliplr(obj.monitorSize);
            scene.position = fliplr(obj.monitorSize)/2;
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            
            p.addStimulus(scene);
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            obj.centerSigmaSelectionIndex = mod(obj.centerSigmaSelectionIndex,length(obj.centerSigmaSequence)) + 1;
            obj.annulusSizeSelectionIndex = mod(obj.annulusSizeSelectionIndex,length(obj.annulusSizeSequence)) + 1;

            obj.rfSigmaCenter = obj.centerSigmaSequence(obj.centerSigmaSelectionIndex);
            obj.rfSigmaSurround = obj.rfSigmaCenter + obj.annulusSizeSequence(obj.annulusSizeSelectionIndex);
        end    
   
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (obj.numberOfAverages*obj.trials);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (obj.numberOfAverages*obj.trials);
        end
        
    end
    
end