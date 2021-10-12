% Spatial split-field spot that switches halves temporally.
classdef splitFieldAdaptation < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Timing
        preTime     = 250 % ms
        stimTime    = 2000 % ms
        tailTime    = 250 % ms
        
        % Stimulus info
        contrasts           = [0 0.9]; % -1 to 1. Light intensities present in split field spot. 
        spotDiameter        = 300; % um
        temporalFrequency   = [60, 30, 20, 15, 12, 10, 6, 5, 3, 2]; % in Hz. How often halves are switched.
        randomize           = true; % randomize between temporal frequencies presented

        % Other details
        rotation            = 0; % deg
        backgroundIntensity = 0.168 % (0-1)
        onlineAnalysis      = 'extracellular'
        numberOfAverages    = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        counter
        order
        timing
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

            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(length(obj.temporalFrequency));
            else
                obj.order = 1:length(obj.temporalFrequency);
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);

            epoch.addParameter('specificTemporalFrequency',obj.temporalFrequency(obj.order(obj.counter+1)));
        end
        
        function p = createPresentation(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            spotDiameterPix = obj.rig.getDevice('Stage').um2pix(obj.spotDiameter);
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity

            %%% Build split-field stimuli
            [xx,yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            th = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2)); % Theta space
            % Adjust for strange monitors
            th = abs(th-pi/2);              
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);
            th = mod(th + obj.rotation + 90,360); % Rotate as needed
            % Define masks
            leftMask = uint8(double(th < 180).*255); % [0 to 180 degrees]
            rightMask = uint8(double(th >= 180).*255); % [180 to 360 degrees]
            %%%%
            
            % Build temporal trajectory
            time = 0:1:obj.stimTime; % in ms
            lightIntensity_Left  = ones(1,length(time)) .* (1+obj.contrasts(1)) .* obj.backgroundIntensity;
            lightIntensity_Right = ones(1,length(time)) .* (1+obj.contrasts(2)) .* obj.backgroundIntensity;
            flipTime        = 0:1000/obj.temporalFrequency(obj.order(obj.counter+1)):obj.stimTime; % in ms
            for a = 2:2:length(flipTime)-1
                range = round(flipTime(a)):round(flipTime(a+1));
                lightIntensity_Left(range)  = (1+obj.contrasts(2)) .* obj.backgroundIntensity;
                lightIntensity_Right(range) = (1+obj.contrasts(1)) .* obj.backgroundIntensity;
            end
            obj.timing = time .* 1e-3; % in s
            
            % Left mask
            left = stage.builtin.stimuli.Rectangle();
            left.position = canvasSize/2;
            left.size = [spotDiameterPix, spotDiameterPix]; 
            leftMaskA = stage.core.Mask(leftMask);
            left.setMask(leftMaskA);
            p.addStimulus(left);
            
            % Define temporal intensity
            leftLightIntensity = stage.builtin.controllers.PropertyController(left,...
                'color', @(state)getLightIntensity(obj, state.time - obj.preTime/1e3, lightIntensity_Left));
            p.addController(leftLightIntensity);
            
            % Right mask
            right = stage.builtin.stimuli.Rectangle();
            right.position = canvasSize/2;
            right.size = [spotDiameterPix, spotDiameterPix]; 
            rightMaskA = stage.core.Mask(rightMask);
            right.setMask(rightMaskA);
            p.addStimulus(right);
            rightLightIntensity = stage.builtin.controllers.PropertyController(right,...
                'color', @(state)getLightIntensity(obj, state.time - obj.preTime .* 1e-3, lightIntensity_Right));
            p.addController(rightLightIntensity);
            
            % Create aperture
            aperture = stage.builtin.stimuli.Rectangle();
            aperture.position = canvasSize/2;
            aperture.color = obj.backgroundIntensity;
            aperture.size = [spotDiameterPix, spotDiameterPix];
            mask = stage.core.Mask.createCircularAperture(1, 1024); %circular aperture
            aperture.setMask(mask);
            p.addStimulus(aperture); %add aperture
            
            % Function for adjusting temporal light intensity
            function s = getLightIntensity(obj, time, lightIntensity)
                if time < 0
                    s = obj.backgroundIntensity;
                elseif time > obj.timing(end)
                    s = lightIntensity(1,end);
                else
                    s = interp1(obj.timing,lightIntensity,time,'nearest');
                end
            end
            
            %hide during pre & post
            grateVisible = stage.builtin.controllers.PropertyController(left, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(grateVisible);
            grateVisible2 = stage.builtin.controllers.PropertyController(right, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(grateVisible2);
            
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages .* length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages .* length(obj.order);
        end
        
    end
    
end