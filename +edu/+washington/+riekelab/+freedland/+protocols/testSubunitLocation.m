% By J. Freedland, 2020.
classdef testSubunitLocation < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime     = 250  % in ms
        stimTime    = 2000 % in ms
        tailTime    = 250  % in ms
        
        % Subunit information
        subunitDiameter         = 40;          % in microns
        subunitLocation_radial  = [25 50 75];  % grid of subunits (in radial coordinates: r, in microns)
        subunitLocation_theta   = [0 120 240]; % grid of subunits (in radial coordinates: theta, in degrees)
        
        % Stimulus parameters
        backgroundIntensity = 0.168;    % intensity of background
        contrast            = 0.9       % contrast
        temporalFrequency   = 4         % Hz
        splitField          = true;
        rotation            = 0;        % degrees
        
        % Additional parameters
        onlineAnalysis      = 'extracellular'
        numberOfAverages    = uint16(3) % number of repeats
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        counter
        order
        masks
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
            
            % Adjust units to DOVES database (units of arcmin)
            subunitRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitDiameter/2,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            subunitLocationPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitLocation_radial,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            videoSize = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                 fliplr(obj.rig.getDevice('Stage').getCanvasSize()),obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            
            % Build subunits
            tempMasks = cell(length(obj.subunitLocation_radial),2);
            for a = 1:length(obj.subunitLocation_radial)
                
                % Build a mask for each subunit
                xCoord = videoSize(2)/2 + subunitLocationPix(a) .* cos(deg2rad(obj.subunitLocation_theta(a)));
                yCoord = videoSize(1)/2 - subunitLocationPix(a) .* sin(deg2rad(obj.subunitLocation_theta(a)));
                [xx,yy] = meshgrid(1:videoSize(2),1:videoSize(1));
                subunit = sqrt((xx - xCoord).^2 + (yy - yCoord).^2) <= subunitRadiusPix;
                
                if obj.splitField == true
                    
                    % Define theta space
                    th = atan((xx - xCoord) ./ (yy - yCoord));
                    th = abs(th-pi/2); 
                    nonsmooth = find(diff(th) > pi/2,1);
                    th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
                    th = rad2deg(th);
                    th = mod(th + obj.rotation,360); % Rotate as required
                    
                    % Split each subunit
                    tempMasks{a,1} = (th > 0 & th <= 180) .* subunit;
                    tempMasks{a,2} = (th > 180 & th <= 360) .* subunit;
                else
                    tempMasks{a,1} = sqrt((xx - xCoord).^2 + (yy - yCoord).^2) <= subunitRadiusPix;
                end
            end
            
            % Turn into full-field masks
            obj.masks{1,1} = zeros(videoSize);
            obj.masks{1,2} = zeros(videoSize);
            for a = 1:size(tempMasks,1)
                obj.masks{1,1} = obj.masks{1,1} + tempMasks{a,1};
                
                if obj.splitField == true
                    obj.masks{1,2} = obj.masks{1,2} + tempMasks{a,2};
                end
            end
            
            obj.counter = 0;
            obj.order = 1:obj.numberOfAverages;
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();

            % Combine masks
            specificDisks = cat(3,obj.masks{1,1},obj.masks{1,2});

            % Add contrast reversing gratings
            for a = 1:2
                grate           = stage.builtin.stimuli.Grating('square'); % Square wave grating
                grate.size      = canvasSize;
                grate.position  = canvasSize / 2;
                grate.spatialFreq = 1 / max(canvasSize * 4); % x2 for diameter, x2 for grating
                grate.color     = 2 * obj.backgroundIntensity; % Amplitude of square wave
                grateShape      = uint8(specificDisks(:,:,a)*255);
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
                spContrast = obj.contrast; % Multiplier on square wave
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
            tf = obj.numEpochsPrepared < length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.order);
        end
    end
end