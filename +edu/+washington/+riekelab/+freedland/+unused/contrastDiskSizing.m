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
        subunitSlices = 4; % slices each radial disk into X quadrants
        
        % Additional options
        randomize = false;
        backgroundIntensity = 0.168 % (0-1)
        rotation = [0 45]; % 0 to 90
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
        radii
        counter
        sequence
        holdData
        specificRotation
        dimTracker
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

            % Convert to pixels
            centerDiskRadii_PIX = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(obj.centerDiskRadii,obj.micronsPerPixel,'UM2PIX'));
            annulusDiskRadii_PIX = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(obj.annulusDiskRadii,obj.micronsPerPixel,'UM2PIX'));
            nearSurroundDiskRadii_PIX = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(obj.nearSurroundDiskRadii,obj.micronsPerPixel,'UM2PIX'));
            farSurroundDiskRadii_PIX = round(max(obj.monitorSize)/2);
            
            % Ensure we have only one vector present
            errorCheck = [length(centerDiskRadii_PIX),length(annulusDiskRadii_PIX),...
                    length(nearSurroundDiskRadii_PIX),length(obj.rotation)];
            if sum(errorCheck > 1) > 1
                error('Please only set one value to be a vector.')
            else
                totalReps = double(max(errorCheck) .* obj.numberOfAverages);
            end
            
            % Collect data
            obj.holdData = [{centerDiskRadii_PIX},{annulusDiskRadii_PIX},...
                    {nearSurroundDiskRadii_PIX},{farSurroundDiskRadii_PIX},{obj.rotation}];    
            
            obj.dimTracker = 1;
            for a = 1:length(obj.holdData)
                if length(obj.holdData{1,a}) == 1
                    obj.holdData{1,a} = repelem(obj.holdData{1,a},totalReps);
                else
                    obj.dimTracker = a;
                    obj.holdData{1,a} = repelem(obj.holdData{1,a},obj.numberOfAverages);
                end
            end
            
            % Account for repeats
            obj.sequence = 1:totalReps;
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
            obj.counter = 1;
        end
        
        function p = createPresentation(obj)
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            
            grabRadii(obj)
            
            % Make disks
            [xx,yy] = meshgrid(1:obj.monitorSize(2),1:obj.monitorSize(1));
            r = sqrt((xx - obj.monitorSize(2)/2).^2 + (yy - obj.monitorSize(1)/2).^2); 
            obj.disks = zeros(obj.monitorSize(1),obj.monitorSize(2),length(obj.radii)-1);
            for a = 1:length(obj.radii)-1
                obj.disks(:,:,a) = r >= obj.radii(a) & r < obj.radii(a+1);
            end
            
            % Add cuts as neccessary
            if obj.subunitSlices > 1
                
                % Define theta space
                th = atan((xx - obj.monitorSize(2)/2) ./ (yy - obj.monitorSize(1)/2));
                th = abs(th-pi/2);              
                nonsmooth = find(diff(th) > pi/2,1);
                th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
                th = rad2deg(th);

                % Region for split field
                sliceAngles = unique(mod((0:360/obj.subunitSlices:360) - obj.specificRotation,360));
                
                overmask = zeros(size(obj.disks,1),size(obj.disks,2),obj.subunitSlices);
                for a = 1:obj.subunitSlices-1
                    overmask(:,:,a) = th >= sliceAngles(a) & th < sliceAngles(a+1);
                end
                overmask(:,:,end) = abs(1 - sum(overmask,3)); % Ignores mod360 annoyances

                % Split masks
                cutDisks = zeros(size(obj.disks,1),size(obj.disks,2),size(overmask,3)*size(obj.disks,3));
                order = 1:size(obj.disks,3);
                counter2 = 1;
                for a = 1:2:obj.subunitSlices
                    for b = order
                        cutDisks(:,:,counter2) = obj.disks(:,:,b) .* overmask(:,:,a);
                        counter2 = counter2 + 1;
                    end
                end
                
                % Switch order so regions are distinct
                altOrder = [order(2:end) order(1)];
                for a = 2:2:obj.subunitSlices
                    for b = altOrder
                        cutDisks(:,:,counter2) = obj.disks(:,:,b) .* overmask(:,:,a);
                        counter2 = counter2 + 1;
                    end
                end
                
                obj.disks = cutDisks;
            end
            
            % Combine disks with same animation (big time saver)
            combinedDisks = zeros(size(obj.disks,1),size(obj.disks,2),2);
            for a = 1:2:size(obj.disks,3)
                combinedDisks(:,:,1) = combinedDisks(:,:,1) + obj.disks(:,:,a);
            end
            for a = 2:2:size(obj.disks,3)
                combinedDisks(:,:,2) = combinedDisks(:,:,2) + obj.disks(:,:,a);
            end

            % Add contrast reversing gratings
            for a = 1:size(combinedDisks,3)
                grate = stage.builtin.stimuli.Grating('square'); % Square wave grating
                grate.size = fliplr(obj.monitorSize);
                grate.position = fliplr(obj.monitorSize)/2;
                grate.spatialFreq = 1 / max(obj.monitorSize * 2);%/(2*obj.diskRadii_PIX(a+1)); % x2 for diameter, x2 for grating
                grate.color = 2*obj.backgroundIntensity; % Amplitude of square wave
                grate.contrast = obj.contrast; % Multiplier on square wave
                
                grateShape = uint8(combinedDisks(:,:,a)*255);
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
            obj.counter = obj.counter + 1;
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            grabRadii(obj);
            epoch.addParameter('specificRotation', obj.specificRotation);
            epoch.addParameter('radii_pixels', obj.radii);
            epoch.addParameter('radii_um', round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(obj.radii,obj.micronsPerPixel,'PIX2UM')));
            epoch.addParameter('tracker',obj.dimTracker);
        end
        
        function grabRadii(obj)
            obj.radii = zeros(1,4);
            for a = 1:4
                obj.radii(a) = obj.holdData{1,a}(obj.counter);
            end
            obj.radii = unique([0 obj.radii]);
            obj.specificRotation = obj.holdData{1,5}(obj.counter);
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.sequence);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.sequence);
        end
        
    end
    
end