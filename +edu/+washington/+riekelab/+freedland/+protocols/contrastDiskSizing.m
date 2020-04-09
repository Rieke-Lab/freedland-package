% Uses cell's RF to find opposed NLs between center/surround
classdef contrastDiskSizing < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms

        % Disk sizing and properties
        centerDiskRadii = [20 30]; % in um.
        annulusDiskRadii = [50 60]; % in um. Set to 0 to ignore.
        nearSurroundDiskRadii = [80 100]; % in um. Set to 0 to ignore.
        contrast = 0.9 % relative to mean (0-1)
        temporalFrequency = 4 % Hz
        
        % Additional options
        randomize = false;
        backgroundIntensity = 0.168 % (0-1)
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(1) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        monitorFrameRate
        micronsPerPixel
        monitorSize
        analysisFigure
        disks
        counter
        radii
        radiiIndex
        correctDim
        sequence
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
            obj.showFigure('edu.washington.riekelab.freedland.figures.contrastDiskSizingFigure',...
                obj.rig.getDevice(obj.amp),'temporalFrequency',obj.temporalFrequency,...
                'preTime',obj.preTime,'stimTime',obj.stimTime);
            
            % Pull variables
            obj.micronsPerPixel = obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.monitorSize = obj.rig.getDevice('Stage').getCanvasSize();
            obj.monitorSize = fliplr(obj.monitorSize); % Adjust to [height, width]

            % Identify disk radii
            centerDiskRadii_PIX = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.centerDiskRadii,obj.micronsPerPixel,'UM2PIX'));
            annulusDiskRadii_PIX = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.annulusDiskRadii,obj.micronsPerPixel,'UM2PIX'));
            nearSurroundDiskRadii_PIX = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.nearSurroundDiskRadii,obj.micronsPerPixel,'UM2PIX'));
            farSurroundDiskRadii_PIX = round(max(obj.monitorSize)/2);
            
            % We can only search along one dimension at a time.
            errorCheck = [length(centerDiskRadii_PIX),length(annulusDiskRadii_PIX),...
                length(nearSurroundDiskRadii_PIX),length(farSurroundDiskRadii_PIX)];
            if sum(errorCheck > 1) > 1
                error('A vector can only be present for the center, annulus, or near-surround radius. Not multiple.')
            end
            
            % Find true radius and ignore regions
            obj.radii = [centerDiskRadii_PIX(1),annulusDiskRadii_PIX(1),nearSurroundDiskRadii_PIX(1),farSurroundDiskRadii_PIX(1)];
            obj.radii(obj.radii == 0) = [];
            obj.radii = [0, obj.radii];
            
            % For sequences to iteratively test
            obj.radiiIndex = find(errorCheck > 1);
            
            if isempty(obj.radiiIndex)
                obj.radiiIndex = 1;
            end

            % Identify correct value to change
            obj.sequence = 1:errorCheck(obj.radiiIndex);
            if obj.radiiIndex <= 1
                obj.correctDim = centerDiskRadii_PIX;
            elseif obj.radiiIndex == 2
                obj.correctDim = annulusDiskRadii_PIX;
            elseif obj.radiiIndex == 3
                obj.correctDim = nearSurroundDiskRadii_PIX;
            end
            obj.radiiIndex = obj.radiiIndex + 1;
            
            % Account for repeats
            obj.sequence = repmat(obj.sequence,1,obj.numberOfAverages);

            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
            obj.counter = 1;
        end
        
        function p = createPresentation(obj)
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            
            % Get next iteration
            specificRadii = obj.correctDim(obj.sequence(obj.counter));
            obj.radii(obj.radiiIndex) = specificRadii;
            obj.counter = obj.counter + 1;
            
            % Make disks
            [xx,yy] = meshgrid(1:obj.monitorSize(2),1:obj.monitorSize(1));
            r = sqrt((xx - obj.monitorSize(2)/2).^2 + (yy - obj.monitorSize(1)/2).^2); 
            obj.disks = zeros(obj.monitorSize(1),obj.monitorSize(2),length(obj.radii)-1);
            for a = 1:length(obj.radii)-1
                obj.disks(:,:,a) = r > obj.radii(a) & r < obj.radii(a+1);
            end

            % Add contrast reversing gratings
            for a = 1:size(obj.disks,3)
                grate = stage.builtin.stimuli.Grating('square'); % Square wave grating
                grate.size = fliplr(obj.monitorSize);
                grate.position = fliplr(obj.monitorSize)/2;
                grate.spatialFreq = 1 / max(obj.monitorSize * 2);%/(2*obj.diskRadii_PIX(a+1)); % x2 for diameter, x2 for grating
                grate.color = 2*obj.backgroundIntensity; % Amplitude of square wave
                grate.contrast = obj.contrast; % Multiplier on square wave
                
                grateShape = uint8(obj.disks(:,:,a)*255);
                grateMask = stage.core.Mask(grateShape);
                grate.setMask(grateMask);
            
                % Contrasting gratings between each radius
                if mod(a,2) == 0
                    grate.phase = 180;
                else
                    grate.phase = 0;
                end

                p.addStimulus(grate); %add grating to the presentation

                grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                    @(state)getGrateContrast(obj, state.time - obj.preTime/1e3));
                p.addController(grateContrast); %add the controller
                
                %hide during pre & post
                grateVisible = stage.builtin.controllers.PropertyController(grate, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(grateVisible);
            end
            
            function c = getGrateContrast(obj, time)
                c = obj.contrast.*sin(2 * pi * obj.temporalFrequency * time);
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            specificRadii = obj.radii;
            changedVal = obj.correctDim(obj.sequence(obj.counter));
            specificRadii(obj.radiiIndex) = changedVal;

            epoch.addParameter('radii_pixels', specificRadii);
            epoch.addParameter('radii_um', round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(specificRadii,obj.micronsPerPixel,'PIX2UM')));
            epoch.addParameter('radii_dimension', obj.radiiIndex);
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.sequence);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.sequence);
        end
        
    end
    
end