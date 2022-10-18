classdef innerCenterFlash < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % ms
        stimTime = 250 % ms
        tailTime = 250 % ms
        innerCenterDiameter = 20:20:300; % inner center, in um
        outerCenterDiameter = 300; % outer center, in um
        gratingsLocation = 'inner center'; % where to place gratings
        gratingsWidth = 20:20:150; % width of gratings, in um
        otherRegionLuminance = 0.319; % light intensity (0-1)
        rotation = 0; % deg
        contrast = 90; % percent contrast of gratings
        backgroundIntensity = 0.168 % (0-1)
        randomizeOrder = true;
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        gratingsLocationType = symphonyui.core.PropertyType('char', 'row', {'inner center', 'outer center', 'none'})
        sequence
        counter
    end
    
    properties (Hidden, Transient)
        analysisFigure
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
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',1);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Arrange experiment conditions
            counter2 = 1;
            
            if strcmp(obj.gratingsLocation,'none')
                obj.sequence = obj.innerCenterDiameter';
                obj.sequence = [obj.sequence, obj.sequence .* 0];
                
                % Remove inner disks > size of outer center
                rmv = (obj.sequence(:,1)) > obj.outerCenterDiameter;
                obj.sequence(rmv,:) = [];
            else
                % Combine inner centers and gratings
                obj.sequence = zeros(length(obj.innerCenterDiameter).*length(obj.gratingsWidth),2);
                for a = 1:length(obj.innerCenterDiameter)
                    for b = 1:length(obj.gratingsWidth)
                        obj.sequence(counter2,:) = [obj.innerCenterDiameter(a), obj.gratingsWidth(b)];
                        counter2 = counter2 + 1;
                    end
                end
                
                if strcmp(obj.gratingsLocation,'inner center')
                    % Remove grating widths > diameter of inner center
                    rmv = (obj.sequence(:,2) .* 2) >= obj.sequence(:,1);
                    obj.sequence(rmv,:) = [];
                elseif strcmp(obj.gratingsLocation,'outer center')
                    % Gratings larger than the outer center
                    rmv = (obj.sequence(:,2) .* 2) >= obj.outerCenterDiameter;
                    obj.sequence(rmv,:) = [];
                end
                
                 % Remove inner disks >= size of outer center 
                 % (= case has no gratings!)
                rmv = (obj.sequence(:,1)) >= obj.outerCenterDiameter; % Can't exactly fill RF
                obj.sequence(rmv,:) = [];
            end
   
            % Randomize order
            if obj.randomizeOrder == true
                obj.sequence = obj.sequence(randperm(size(obj.sequence,1)),:);
            end

            obj.counter = 1;
            obj.sequence = repmat(obj.sequence,obj.numberOfAverages,1);
            
            %%% Quick bug fix: very last epoch won't play, so we double last condition
            obj.sequence = [obj.sequence; obj.sequence(end,:)];
        end
             
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('currentInnerCenterDiameter', obj.sequence(obj.counter,1));
            epoch.addParameter('currentGratingsWidth', obj.sequence(obj.counter,2));
        end

        function p = createPresentation(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            
            % Convert from microns to pixels
            innerDiameterPix = obj.rig.getDevice('Stage').um2pix(obj.sequence(obj.counter,1));
            gratingWidthPix   = obj.rig.getDevice('Stage').um2pix(obj.sequence(obj.counter,2));
            centerDiameterPix = obj.rig.getDevice('Stage').um2pix(obj.outerCenterDiameter);
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            
            % Inner center mask
            inner  = stage.builtin.stimuli.Grating('square'); %square wave grating
            inner.orientation = obj.rotation;
            inner.position    = canvasSize/2;
            inner.color       = obj.backgroundIntensity .* 2;
            inner.contrast    = obj.contrast/100;
            inner.size = [innerDiameterPix, innerDiameterPix];
            inner.spatialFreq = 1/(2*gratingWidthPix); %convert from bar width to spatial freq
            mask = stage.core.Mask.createAnnulus(0, 1, 1024); %circular aperture
            inner.setMask(mask);
            
            % Outer center mask
            outer  = stage.builtin.stimuli.Grating('square'); %square wave grating
            outer.orientation = obj.rotation;
            outer.position    = canvasSize/2;
            outer.color       = obj.backgroundIntensity .* 2;
            outer.contrast    = obj.contrast/100;
            outer.size = [centerDiameterPix, centerDiameterPix];
            outer.spatialFreq = 1/(2*gratingWidthPix); %convert from bar width to spatial freq
            mask = stage.core.Mask.createAnnulus(innerDiameterPix/centerDiameterPix, 1, 1024); %circular aperture
            outer.setMask(mask);
            
            % Only apply gratings where specified
            if strcmp(obj.gratingsLocation,'inner center')
                
                % Outer center becomes uniform disk
                outer.spatialFreq = 0;
                outer.color = obj.otherRegionLuminance .* 2;
                p.addStimulus(outer);
                
                % Build contrast reversing boundary
                zeroCrossings   = 0:(inner.spatialFreq^-1):inner.size(1); 
                offsets         = zeroCrossings-inner.size(1)/2; %difference between each zero crossing and center of texture, pixels
                [shiftPix, ~]   = min(offsets(offsets>0)); %positive shift in pixels
                phaseShift_rad  = (shiftPix/(inner.spatialFreq^-1))*(2*pi); %phaseshift in radians
                phaseShift      = 360*(phaseShift_rad)/(2*pi); %phaseshift in degrees
                inner.phase = phaseShift; 
                p.addStimulus(inner);
                
            elseif strcmp(obj.gratingsLocation,'outer center')
                
                % Inner center becomes uniform disk
                inner.spatialFreq = 0;
                inner.color = obj.otherRegionLuminance .* 2;
                p.addStimulus(inner);
               
                % Build contrast reversing boundary
                zeroCrossings   = 0:(outer.spatialFreq^-1):outer.size(1); 
                offsets         = zeroCrossings-outer.size(1)/2; %difference between each zero crossing and center of texture, pixels
                [shiftPix, ~]   = min(offsets(offsets>0)); %positive shift in pixels
                phaseShift_rad  = (shiftPix/(outer.spatialFreq^-1))*(2*pi); %phaseshift in radians
                phaseShift      = 360*(phaseShift_rad)/(2*pi); %phaseshift in degrees
                outer.phase     = phaseShift;
                p.addStimulus(outer);
                
            elseif strcmp(obj.gratingsLocation,'none')
                inner.spatialFreq = 0;
                outer.spatialFreq = 0;
                
                inner.color = obj.otherRegionLuminance .* 2;
                outer.color = obj.otherRegionLuminance .* 2;
                
                p.addStimulus(inner);
                p.addStimulus(outer);
            end
                
            innerVisible = stage.builtin.controllers.PropertyController(inner, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(innerVisible);
            
            outerVisible = stage.builtin.controllers.PropertyController(outer, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(outerVisible);
            
            obj.counter = obj.counter + 1;
        end
 
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < (size(obj.sequence,1)-1); % Present one less epoch
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < (size(obj.sequence,1)-1);
        end
        
        
    end
    
end