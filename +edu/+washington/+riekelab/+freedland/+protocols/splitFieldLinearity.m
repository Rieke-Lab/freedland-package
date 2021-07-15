% Tests whether post-rectification split-field flashes perform similarly to
% uniform disks

classdef splitFieldLinearity < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Timing
        preTime     = 250 % ms
        stimTime    = 250 % ms
        tailTime    = 250 % ms
        
        % Stimulus info
        leftContrasts   = 0:0.2:1 % on left side of split gratings (-1 to 1)
        rightContrasts  = [-0.75 -0.25]% on right side of split gratings (-1 to 1)
        spotDiameter    = 300; % um
        
        % Control
        rectificationThreshold  = 0; % rectify values beyond contrast value
        cellClass               = 'ON' % type of cell
        includeUniformFlashes   = true; % include uniform flashes with expected post-rectification luminance

        % Other details
        rotation = 0; % deg
        backgroundIntensity = 0.168 % (0-1)
        randomize = true; % random order of presentation
        onlineAnalysis = 'none'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        cellClassType = symphonyui.core.PropertyType('char', 'row', {'ON', 'OFF'})   
        leftLuminances
        rightLuminances
        counter
        order
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
            
            % All possible combinations
            obj.leftLuminances = repelem(obj.leftContrasts,length(obj.rightContrasts));
            obj.rightLuminances = repmat(obj.rightContrasts,1,length(obj.leftContrasts));
            
            if obj.includeUniformFlashes == true
                L = obj.leftLuminances;
                R = obj.rightLuminances;
                if strcmp(obj.cellClass,'ON')
                    L(L < obj.rectificationThreshold) = obj.rectificationThreshold;
                    R(R < obj.rectificationThreshold) = obj.rectificationThreshold;
                else
                    L(L > obj.rectificationThreshold) = obj.rectificationThreshold;
                    R(R > obj.rectificationThreshold) = obj.rectificationThreshold;
                end
                integratedLuminances = nanmean([L;R],1);
                
                obj.leftLuminances = [obj.leftLuminances, integratedLuminances];
                obj.rightLuminances = [obj.rightLuminances, integratedLuminances];
            end
            
            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(length(obj.leftLuminances));
            else
                obj.order = 1:length(obj.leftLuminances);
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('specificLeftContrast',obj.leftContrasts(obj.order(obj.counter+1)));
            epoch.addParameter('specificRightContrast',obj.rightContrasts(obj.order(obj.counter+1)));
        end
        
        function p = createPresentation(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            
            % Convert from microns to pixels
            spotDiameterPix = obj.rig.getDevice('Stage').um2pix(obj.spotDiameter);
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            
            % Define masks
            [xx,yy] = meshgrid(1:canvasSize(2),1:canvasSize(1));
            th = atan((xx - canvasSize(2)/2) ./ (yy - canvasSize(1)/2));
            th = abs(th-pi/2);              
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);
            th = mod(th + obj.rotate,360); % Rotate as required
            leftMask = uint8(double(th < 180) .* 255);
            rightMask = uint8(double(th >= 180) .* 255);
            
            % Convert contrast to luminance
            leftL = (obj.leftLuminances(obj.order(obj.counter+1)) + 1).*obj.backgroundIntensity;
            rightL = (obj.rightLuminances(obj.order(obj.counter+1)) + 1).*obj.backgroundIntensity;
            
            % Left mask
            left = stage.builtin.stimuli.Grating();
            left.color = 2*leftL; % x2 for grating preset
            left.position = canvasSize/2;
            left.size = [spotDiameterPix, spotDiameterPix]; 
            leftMaskA = stage.core.Mask(leftMask);
            left.setMask(leftMaskA);
            
            % Right mask
            right = stage.builtin.stimuli.Grating();
            right.color = 2*rightL; % x2 for grating preset
            right.position = canvasSize/2;
            right.size = [spotDiameterPix, spotDiameterPix]; 
            rightMaskA = stage.core.Mask(rightMask);
            left.setMask(rightMaskA);
            
            % Create aperture
            aperture = stage.builtin.stimuli.Rectangle();
            aperture.position = canvasSize/2;
            aperture.color = obj.backgroundIntensity;
            aperture.size = [spotDiameterPix, spotDiameterPix];
            mask = stage.core.Mask.createCircularAperture(1, 1024); %circular aperture
            aperture.setMask(mask);
            p.addStimulus(aperture); %add aperture
            
            %hide during pre & post
            grateVisible = stage.builtin.controllers.PropertyController(left, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(grateVisible);
            
            grateVisible2 = stage.builtin.controllers.PropertyController(right, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(grateVisible2);
            
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end
        
        %same presentation each epoch in a run. Replay.
        function controllerDidStartHardware(obj)
            controllerDidStartHardware@edu.washington.riekelab.protocols.RiekeLabProtocol(obj);
            if (obj.numEpochsCompleted >= 1) && (obj.numEpochsCompleted < obj.numberOfAverages)
                obj.rig.getDevice('Stage').replay
            else
                obj.rig.getDevice('Stage').play(obj.createPresentation());
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