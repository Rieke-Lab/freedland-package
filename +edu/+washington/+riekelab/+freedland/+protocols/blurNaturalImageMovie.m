% Plays a natural image movie, courtesy of the DOVES database.
classdef blurNaturalImageMovie < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        
        % Stimulus timing
        preTime     = 250   % in ms
        stimTime    = 5500  % in ms
        tailTime    = 250   % in ms
        
        % Natural image trajectory
        imageNo     = 71;    % natural image number (1 to 101)
        observerNo  = 1;    % observer number (1 to 19)
        
        % Blur information: how to blur each sequential step?
        includeNaturalMovie = true; % whether to include unblurred variant
        coneBlur    = 1.5; % sigma of Gaussian blur kernel (in microns) - first stage of filtering
        subunitBlur = 15; % sigma of Gaussian blur kernel (in microns) - second stage of filtering
        lowerRectification = [-30, -15, 0, 15, Inf]; % Rectify values below each value. Inf ignores rectification.
        upperRectification = [-30, -15, 0, 15, Inf];  % Rectify values above each value. Inf ignores rectification.
        rgcBlur	= [0 50 75 100]; % sigma of Gaussian blur kernel (in microns) - last stage of filtering

        % Set region for testing
        centerDiameter = 300; % only present natural image in RF center (in microns). Set to 0 to ignore.

        % Set image type
        naturalImage = 'natural'; % natural image or select Fourier statistics.
        
        % Additional parameters
        randomize        = true; % whether to randomize movies shown
        onlineAnalysis   = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
        
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        naturalImageType = symphonyui.core.PropertyType('char', 'row', {'natural', 'phase', 'magnitude'});
        backgroundIntensity
        xTraj
        yTraj
        timeTraj
        imageMatrix
        sequence
        counter
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
            
            % Gather natural image information
            [path,image] = edu.washington.riekelab.freedland.videoGeneration.utils.pathDOVES(obj.imageNo, obj.observerNo);
            image = image ./ max(image(:)); % Normalize (scale to monitor)
            obj.backgroundIntensity = mean(image(:));
            obj.imageMatrix = uint8(image .* 255);
                    
            % Isolate DOVES eye trajectories
            obj.xTraj = path.x;
            obj.yTraj = path.y;
            
            % Invert for monitor and adjust position relative to center of image
            obj.xTraj = -(obj.xTraj - size(image,2)/2);
            obj.yTraj = (obj.yTraj - size(image,1)/2);
            
            % Convert from DOVES units (1 px = 1 arcmin) to monitor units
            obj.xTraj = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.xTraj,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
            obj.yTraj = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.yTraj,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix');
            obj.timeTraj = (0:(length(obj.xTraj)-1)) ./ 200; % convert DOVES resolution (200Hz) to seconds
            
            % Define all combinations manually (as a tree)
            obj.sequence = [];
            for a = 1:length(obj.coneBlur)
                for b = 1:length(obj.subunitBlur)
                    for c = 1:length(obj.lowerRectification)
                        for d = 1:length(obj.upperRectification)
                            for f = 1:length(obj.rgcBlur)
                                tmp = [obj.coneBlur(a), obj.subunitBlur(b), obj.lowerRectification(c),...
                                       obj.upperRectification(d) obj.rgcBlur(f)];
                                obj.sequence = [obj.sequence; tmp];
                            end
                        end
                    end
                end
            end
            
            if obj.includeNaturalMovie == true
                % +Inf == raw movie condition
                obj.sequence = [obj.sequence; repelem(Inf, 1, size(obj.sequence,2))];
            end
            
            if strcmp(obj.naturalImage,'phase') || strcmp(obj.naturalImage,'magnitude')
                % -Inf == raw movie condition (with only isolated Fourier statistic)
                obj.sequence = [obj.sequence; repelem(-Inf, 1, size(obj.sequence,2))]; % Raw Fourier condition
            end
            
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(size(obj.sequence,1)),:);
            end
            obj.counter = 0;
            
            t = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3 + 1.5; % time (+1.5s rig delay)
            disp(strcat('est. protocol time:',mat2str(length(obj.sequence) .* obj.numberOfAverages .* t / 60),'min'));
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;

            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            % Export specific parameters
            epoch.addParameter('coneBlur_specific',obj.sequence(obj.counter+1,1));
            epoch.addParameter('subunitBlur_specific',obj.sequence(obj.counter+1,2));
            epoch.addParameter('lowerRectification_specific',obj.sequence(obj.counter+1,3));
            epoch.addParameter('upperRectification_specific',obj.sequence(obj.counter+1,4));
            epoch.addParameter('rgcBlur_specific',obj.sequence(obj.counter+1,5));
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize(); % in normal pixels            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity)

            %%% Filter image
            % Units are in arcmin - we filter in the image at the scale of
            % the DOVES database, then enlarge to fit the monitor
            selection = obj.sequence(obj.counter+1,:);
            disp(selection)
            coneBlur_arcmin = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(selection(1),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            subunitBlur_arcmin = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(selection(2),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            rgcBlur_arcmin = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(selection(5),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            
            % Pull raw image
            tmp = double(obj.imageMatrix);
            
            % Isolate Fourier properties
            if strcmp(obj.naturalImage,'phase') || strcmp(obj.naturalImage,'magnitude')
                if sum(selection == Inf) ~= length(selection)
                    
                    % Isolate amplitude, phase
                    FFT = fftshift(fft2(tmp));
                    amplitude = abs(FFT);
                    phase = unwrap(angle(FFT));
                
                    % Sanity check
%                     originalImage = abs(ifft2(amplitude .* cos(phase) + amplitude .* sin(phase) .* 1i));

                    if strcmp(obj.naturalImage,'phase')
                        tmp = abs(ifft2(cos(phase) + sin(phase) .* 1i));

                        % Adjust background intensity to match original image
                        for iter = 1:10
                            tmp = tmp ./ max(tmp(:)) .* 255;
                            tmp = tmp + (obj.backgroundIntensity - nanmean(tmp(:))); 
                        end
                    elseif strcmp(obj.naturalImage,'magnitude')
                        tmp = abs(ifft2(amplitude));
                        tmp(tmp > 255) = 255;
                        
                        % Adjust background intensity to match original image
                        tmp = tmp + (obj.backgroundIntensity - nanmean(tmp(:)));
                    end
                end
            end
            
            % Apply blur
            if sum(selection == Inf) ~= length(selection) && ... % All inf = unblurred condition
               sum(selection == -Inf) ~= length(selection)       % All -inf = unblurred Fourier condition

                % Apply cone blur
                tmp = imgaussfilt(tmp,coneBlur_arcmin);
                
                % Apply subunit blur
                tmp = imgaussfilt(tmp,subunitBlur_arcmin);
                
                %%% Rectify
                % Convert rectification bounds from % contrast to 8-bit luminance
                l_limit = (selection(3)/100 + 1) .* obj.backgroundIntensity .* 255;
                u_limit = (selection(4)/100 + 1) .* obj.backgroundIntensity .* 255;
                if sum(selection(3:4) == Inf) == 0 % Rectification = inf = nonrectified condition
                    tmp(tmp < l_limit) = l_limit;
                    tmp(tmp > u_limit) = u_limit;
                end
                
                % Apply RGC blur
                if rgcBlur_arcmin > 0
                    tmp = imgaussfilt(tmp,rgcBlur_arcmin);
                end
            end
            
            % Insert image and sizing information for stage.
            scene = stage.builtin.stimuli.Image(uint8(tmp));
            scene.size = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(fliplr(size(obj.imageMatrix)),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'arcmin2pix'); % Convert to monitor units
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            % Control the position of the image as a function of time.
            scenePosition = stage.builtin.controllers.PropertyController(scene,...
                'position', @(state)getScenePosition(obj, state.time - obj.preTime/1e3, p0));

            function p = getScenePosition(obj, time, p0)
                if time <= 0 % Before stimulus begins
                    p = p0;
                elseif time > obj.timeTraj(end) % Beyond eye trajectory, hang on last frame
                    p(1) = p0(1) + obj.xTraj(end);
                    p(2) = p0(2) + obj.yTraj(end);
                else % Within eye trajectory and stim time
                    dx = interp1(obj.timeTraj,obj.xTraj,time);
                    dy = interp1(obj.timeTraj,obj.yTraj,time);
                    p(1) = p0(1) + dx;
                    p(2) = p0(2) + dy;
                end
            end
            
            % Add information to Stage
            p.addStimulus(scene);
            p.addController(scenePosition);
            
            % Add additional controller: mediate when stimulus is visible.
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            % Create center mask
            if (obj.centerDiameter > 0)
                
                % Convert to monitor pixels
                maskDiameterPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(obj.centerDiameter,...
                    obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
                
                aperture = stage.builtin.stimuli.Rectangle();
                aperture.position = canvasSize/2;
                aperture.color = obj.backgroundIntensity;
                aperture.size = [max(canvasSize) max(canvasSize)];
                mask = stage.core.Mask.createCircularAperture(maskDiameterPix/max(canvasSize), 1024); %circular aperture
                aperture.setMask(mask);
                p.addStimulus(aperture); %add aperture
            end
            
            obj.counter = mod(obj.counter + 1,size(obj.sequence,1));
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < size(obj.sequence,1) .* obj.numberOfAverages;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < size(obj.sequence,1) .* obj.numberOfAverages;
        end
    end
end
