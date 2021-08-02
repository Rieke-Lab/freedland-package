% Temporally tests correlation structures in low-dimensional subunit space
classdef subunitCorrelation < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime     = 250;   % ms
        stimTime    = 2000;  % ms
        tailTime    = 250;   % ms

        % Disk sizing and properties
        centerDiameter      = 200   % in um
        centerDivisions     = 8;    % number of spatial divisions in RF center
        rotate              = 0;    % rotate spatial divisions (in deg)
        contrast            = 0.9;  % temporal changes in luminance (0 to 1)
        backgroundIntensity = 0.168; % light intensity of background (0 to 1)
        temporalFrequency   = 4;    % temporal behavior (Hz)
        
        % Additional options
        randomize           = true; % randomize order of epochs
        onlineAnalysis      = 'extracellular'
        numberOfAverages    = uint16(5) % number of repeats to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        masks
        positiveContrast
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

            % Convert units from microns to pixels
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            centerRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.centerDiameter/2,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            
            % Create radial space
            [xx,yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            r = sqrt((xx - canvasSize(1)/2).^2 + (yy - canvasSize(2)/2).^2); % Rho space
            th = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2)); % Theta space
            th = abs(th-pi/2);              
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);
            th = mod(th + obj.rotate,360); % Rotate as required
            if obj.centerDivisions == 0
                theta = [0,360];
            else
                theta = 0:360/obj.centerDivisions:360;
            end
            
            % Generate spatial masks
            obj.masks = zeros([fliplr(canvasSize),length(theta)-1]);
            for a = 1:length(theta) - 1
                obj.masks(:,:,a) = (th >= theta(a) & th < theta(a+1)) .* (r <= centerRadiusPix);
            end
            
            % Masks that contain positive stimulus
            obj.positiveContrast = nchoosek(1:size(obj.masks,3),ceil(size(obj.masks,3)/2));
            
            expectedTime = (obj.preTime + obj.stimTime + obj.tailTime + 1.5*1000) / (1000*60); % +1.5 sec rig delay, in min
            disp(strcat('Approx. stimulus time: ',mat2str(round(size(obj.positiveContrast,1)*expectedTime*obj.numberOfAverages)),'min'));

            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(size(obj.positiveContrast,1));
            else
                obj.order = 1:size(obj.positiveContrast,1);
            end
        end
        
        function prepareEpoch(obj, epoch)

            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('alignedPositiveContrastRegions',obj.positiveContrast(obj.order(obj.counter+1),:)); % About center (for figure)
        end
        
        function p = createPresentation(obj)

            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();

            regions = obj.positiveContrast(obj.order(obj.counter+1),:);
            ordering = false(1,size(obj.masks,3));
            ordering(regions) = true;
            summedMasks(:,:,1) = sum(obj.masks(:,:,ordering),3);
            summedMasks(:,:,2) = sum(obj.masks(:,:,~ordering),3);

            % Add contrast reversing gratings
            for a = 1:2
                grate           = stage.builtin.stimuli.Grating('square'); % Square wave grating
                grate.size      = canvasSize;
                grate.position  = canvasSize / 2;
                grate.spatialFreq = 1 / max(canvasSize * 4); % x2 for diameter, x2 for grating
                grate.color     = 2 * obj.backgroundIntensity; % Amplitude of square wave
                spContrast      = obj.dotContrast; % Multiplier on square wave
                grateShape      = uint8(summedMasks(:,:,a)*255);
                grateMask       = stage.core.Mask(grateShape);
                grate.setMask(grateMask);
            
                % Contrasting gratings between each set of masks
                if a == 1
                    grate.phase = 180;
                else
                    grate.phase = 0;
                end

                p.addStimulus(grate); % Add grating to the presentation

                % Control contrast
                grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                    @(state)getGrateContrast(obj, spContrast, state.time - obj.preTime/1e3));
                p.addController(grateContrast); % Add the controller
                
                % Hide during pre & post
                grateVisible = stage.builtin.controllers.PropertyController(grate, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(grateVisible);
            end
            obj.counter = mod(obj.counter + 1,length(obj.order));

            function c = getGrateContrast(obj, specificContrast, time)
                c = specificContrast.*sin(2 * pi * obj.temporalFrequency * time);
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.order) .* obj.numberOfAverages;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.order) .* obj.numberOfAverages;
        end 
        
    end
    
end