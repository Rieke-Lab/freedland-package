% Performs contrast-reversing stimuli in polar space
classdef contrastReversingSpatialCenter < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        % Basic parameters
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms

        % Receptive field settings
        sigmaCenter = 60; % center gaussian from receptive field fit
        sigmaSurround = 170;% surround gaussian from receptive field fit

        % Spatial parameters
        wedges = 8; % radial wedges
        innerCenter = 40; % percent of inner RF
        outerCircles = 3; % number of circles to impose in outer RF

        % Stimulus parameters
        temporalFrequency = 4 % Hz
        backgroundIntensity = 0.168 % (0-1)

        % Additional
        randomizeOrder = true;
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        masks
        maskTracker
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
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));

            %%% Build associated masks
            % Convert from microns to pixels
            centerRF_px     = obj.rig.getDevice('Stage').um2pix(obj.sigmaCenter);
            surroundRF_px   = obj.rig.getDevice('Stage').um2pix(obj.sigmaSurround);

            % Generate 2D Gaussians
            canvasSize      = obj.rig.getDevice('Stage').getCanvasSize();
            centerGauss     = fspecial('gaussian',fliplr(canvasSize),centerRF_px);
            surroundGauss   = fspecial('gaussian',fliplr(canvasSize),surroundRF_px);
        
            % Calculate difference of gaussians (receptive field)
            diffGaussian    = centerGauss - surroundGauss;
            RFFilter        = diffGaussian ./ max(diffGaussian(:)); % Normalize filter

            % Isolate receptive field center
            RFFilter    = RFFilter .* (RFFilter > 0);

            % Define polar space
            [xx,yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            r = sqrt((xx - canvasSize(1)/2).^2 + (yy - canvasSize(2)/2).^2);
            th = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2));
            th = abs(th-pi/2);              
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);

            % Integrate radially
            radii = 0:1:max(r(RFFilter > 0));
            cumulativeSum = zeros(length(radii),1);
            for a = 1:length(radii)
                cumulativeSum(a) = sum(RFFilter .* (r < radii(a)),'all');
            end

            % Isolate inner center (via total integrated area)
            innerCutoff = max(cumulativeSum) .* obj.innerCenter/100;
            [~,i]       = min(abs(cumulativeSum - innerCutoff));
            innerMask   = (r < radii(i));

            % Identify outer center
            outerMask   = zeros([size(innerMask),obj.outerCircles]);
            resolution  = (100 - obj.innerCenter) ./ obj.outerCircles;
            for a = 1:obj.outerCircles
                cutoff  = max(cumulativeSum) .* (obj.innerCenter + resolution .* a)/100;
                [~,i]   = min(abs(cumulativeSum - cutoff));
                outerMask(:,:,a) = (r < radii(i)) .* (RFFilter > 0); % Isolate center

                % Remove other masks
                outerMask(:,:,a) = outerMask(:,:,a) .* (innerMask == 0) .* (sum(outerMask,3) == 1);
            end

%             % Sanity check: confirm areas are appropriately-sized
%             A = [squeeze(sum(innerMask .* RFFilter,[1 2])); ...
%                  squeeze(sum(outerMask .* RFFilter,[1 2]))];
%             A = A ./ max(cumulativeSum) .* 100; % Percent of integrated RF contained in each mask

            % Turn all circular masks into wedges
            obj.wedges(obj.wedges == 0) = 1;
            theta = 0:360/obj.wedges:360;
            tmp = cat(3,innerMask,outerMask);
            m = cell(obj.wedges,size(tmp,3));
            for b = 1:size(tmp,3)
                for a = 1:length(theta) - 1
                    radialMask = th >= theta(a) & th < theta(a+1); % Angular filter (theta)
                    m{a,b} = radialMask .* tmp(:,:,b); % Individual masks
                end
            end

            %%% Build unique sets of contrast-reversing stimuli
            obj.masks = [];
            obj.maskTracker = [];
            
            wedgeMasks = [];
            counter = 1;

            % Build cases with only wedges
            for condition = 1:2 % [uniform, checkerboard]
                % Build stimulus
                tmp = zeros(fliplr(canvasSize));
                for a = 1:size(m,2)
                    for b = [1 3 5 7]
                        % Apply checkerboard
                        if condition == 2 & mod(a,2) == 1
                            b = b+1;
                        end
                        tmp = tmp + m{b,a};
                    end
                end
                wedgeMasks{condition,1} = (tmp == 1);
            end

            %%% Cases with only circles
            for k = 0:obj.outerCircles
                tmp = nchoosek(1:obj.outerCircles,k);
                for a = 1:size(tmp,1)
                    t = zeros(1,obj.outerCircles);
                    t(tmp(a,:)) = 1;

                    % Define mask
                    onset = sum(outerMask(:,:,t == 1),3);

                    if sum(t,2) > 1
                        t2 = 1:size(wedgeMasks,1);
                    else
                        t2 = 1;
                    end

                    % Combine with wedge masks
                    for b = t2
                        obj.masks{counter,1} = wedgeMasks{b,1} .* onset;
                        obj.masks{counter,2} = (obj.masks{counter,1} ~= 1) .* (RFFilter > 0);
                        counter = counter+1;

                        % Include condition with spatial structure in center
                        obj.masks{counter,1} = obj.masks{counter-1,1} + wedgeMasks{b,1} .* innerMask;
                        obj.masks{counter,2} = (obj.masks{counter,1} ~= 1) .* (RFFilter > 0);
                        counter = counter+1;
                    end
                end
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
        end

        function p = createPresentation(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            
            %convert from microns to pixels...
            apertureDiameterPix = obj.rig.getDevice('Stage').um2pix(obj.apertureDiameter);
            centerBarWidthPix   = obj.rig.getDevice('Stage').um2pix(obj.sequence(obj.counter,2));
            surroundBarWidthPix = obj.rig.getDevice('Stage').um2pix(obj.sequence(obj.counter,4));
            centerContrast_sp = obj.sequence(obj.counter,1);
            surroundContrast_sp = obj.sequence(obj.counter,3);
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            
            % Create grating in RF center and surround
            for a = 1:2 % Center first, then surround
                grate  = stage.builtin.stimuli.Grating('square'); %square wave grating
                grate.orientation = obj.rotation;
                grate.position    = canvasSize/2;
                grate.color       = 2*obj.backgroundIntensity;
                
                if a == 1 % Center
                    grate.size = [apertureDiameterPix, apertureDiameterPix];
                    grate.spatialFreq = 1/(2*centerBarWidthPix); %convert from bar width to spatial freq
                elseif a == 2 % Surround
                    grate.size = [max(canvasSize) max(canvasSize)];
                    if surroundBarWidthPix > 0
                        grate.spatialFreq = 1/(2*surroundBarWidthPix); %convert from bar width to spatial freq
                    end
                    % Apply circular hole for center
                    mask = stage.core.Mask.createCircularAperture(apertureDiameterPix/max(canvasSize), 1024); %circular aperture
                    grate.setMask(mask);
                end

                %calc to apply phase shift s.t. a contrast-reversing boundary
                %is in the center regardless of spatial frequency. Arbitrarily
                %say boundary should be positve to right and negative to left
                %crosses x axis from neg to pos every period from 0
                zeroCrossings   = 0:(grate.spatialFreq^-1):grate.size(1); 
                offsets         = zeroCrossings-grate.size(1)/2; %difference between each zero crossing and center of texture, pixels
                [shiftPix, ~]   = min(offsets(offsets>0)); %positive shift in pixels
                phaseShift_rad  = (shiftPix/(grate.spatialFreq^-1))*(2*pi); %phaseshift in radians
                phaseShift      = 360*(phaseShift_rad)/(2*pi); %phaseshift in degrees
                grate.phase = phaseShift; 
                p.addStimulus(grate);

                % Make it contrast-reversing
                if (obj.temporalFrequency > 0) 
                    if a == 1 % Center
                        grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                            @(state)getGrateContrastCenter(obj, state.time - obj.preTime/1e3, centerContrast_sp));
                    else
                        grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                            @(state)getGrateContrastSurround(obj, state.time - obj.preTime/1e3, surroundContrast_sp));
                    end
                    p.addController(grateContrast); %add the controller
                end
                
                grateVisible = stage.builtin.controllers.PropertyController(grate, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(grateVisible);
            end
            
            function c = getGrateContrastCenter(obj, time, contrast)
                c = contrast.*sin(2 * pi * obj.temporalFrequency * time);
            end
            
            function c = getGrateContrastSurround(obj, time, contrast)
                if strcmp(obj.surroundPhase,'in-phase')
                    c = contrast.*sin(2 * pi * obj.temporalFrequency * time); % -1 to counter (plays after adding +1)
                elseif strcmp(obj.surroundPhase,'out-of-phase')
                    c = contrast.*sin(2 * pi * obj.temporalFrequency * time) .* -1;
                end
            end
            
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