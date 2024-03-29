% Plays a natural image movie, courtesy of the DOVES database.
classdef blurNaturalImageMovie < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        
        % Stimulus timing
        preTime     = 250   % in ms
        stimTime    = 5500  % in ms
        tailTime    = 250   % in ms
        
        % Natural image trajectory
        imageNo     = 5;    % natural image number (1 to 101)
        observerNo  = 1;    % observer number (1 to 19)
        
        % Blur information: how to blur each sequential step?
        includeNaturalMovie = true; % whether to include unblurred variant
        coneBlur    = 1.5; % sigma of Gaussian blur kernel (in microns) - first stage of filtering
        subunitBlur = 15; % sigma of Gaussian blur kernel (in microns) - second stage of filtering
        lowerRectification = [-90, -60, -30, 0, 30, Inf]; % Rectify values below each value. Inf ignores rectification.
        upperRectification = [-30, 0, 30, 60, 90, Inf]; % Rectify values above each value. Inf ignores rectification.
        rectificationSite  = 'post-subunit'; % where to place rectifier (both = test pre- and post-subunit rectification individually).
        rgcBlur     = 0; % sigma of Gaussian blur kernel (in microns) - last stage of filtering
        
        % Binning pixels
        binPixels = true; % bins pixels into values.
        binSize = [0 10 20 30 50]; % in % contrast

        % Set region for testing
        centerDiameter = 300; % only present natural image in RF center (in microns). Set to 0 to ignore.

        % Additional parameters
        randomize        = true; % whether to randomize movies shown
        onlineAnalysis   = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
        
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        rectificationSiteType = symphonyui.core.PropertyType('char', 'row', {'pre-subunit', 'post-subunit', 'testBoth'});
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

            % Identify where to place rectification threshold
            if strcmp(obj.rectificationSite,'pre-subunit') % 0
                obj.sequence = [obj.sequence,zeros(size(obj.sequence,1),1)];
            elseif strcmp(obj.rectificationSite,'post-subunit') % 1
                obj.sequence = [obj.sequence,ones(size(obj.sequence,1),1)];
            elseif strcmp(obj.rectificationSite,'testBoth')
                obj.sequence = [obj.sequence,zeros(size(obj.sequence,1),1)];
                obj.sequence = [obj.sequence;[obj.sequence(:,1:end-1), ones(size(obj.sequence,1),1)]];
            end
            
            % Set up binning conditions
            if obj.binPixels == false
                obj.sequence = [obj.sequence zeros(size(obj.sequence,1),1)];
            else
                r = repelem(obj.binSize',size(obj.sequence,1),1);
                obj.sequence = [repmat(obj.sequence,length(obj.binSize),1) r];
            end
            
            if obj.includeNaturalMovie == true
                % +Inf == raw movie condition
                obj.sequence = [repelem(Inf, 1, size(obj.sequence,2)); obj.sequence];
            end
            
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(size(obj.sequence,1)),:);
            end
            obj.counter = 0;
            
            t = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3 + 1.5; % time (+1.5s rig delay)
            disp(strcat('est. protocol time:',mat2str(length(obj.sequence) .* obj.numberOfAverages .* t / 60),'min'));
            disp('...')
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;

            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            % Export specific parameters
            A = obj.sequence(obj.counter+1,:);
            epoch.addParameter('coneBlur_specific',A(1));
            epoch.addParameter('subunitBlur_specific',A(2));
            epoch.addParameter('lowerRectification_specific',A(3));
            epoch.addParameter('upperRectification_specific',A(4));
            epoch.addParameter('rgcBlur_specific',A(5));
            
            rectificationSites = {'pre-subunit', 'post-subunit'};
            if A(6) ~= Inf % Non-natural image condition
                epoch.addParameter('rectificationSite_specific',rectificationSites{A(6)+1});
            else
                epoch.addParameter('rectificationSite_specific',Inf);
            end
            epoch.addParameter('pixelBinning_specific',A(7));
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
            coneBlur_arcmin = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(selection(1),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            subunitBlur_arcmin = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(selection(2),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            rgcBlur_arcmin = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(selection(5),...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2arcmin');
            
            % Pull raw image
            tmp = double(obj.imageMatrix);
            
            % Apply blur
            if sum(selection == Inf) ~= length(selection) % All inf = natural image condition

                % Apply cone blur
                if coneBlur_arcmin > 0
                    tmp = imgaussfilt(tmp,coneBlur_arcmin);
                    disp(strcat(mat2str(selection(1)),'um cone blur applied...'))
                end

                % Convert rectification bounds from % contrast to 8-bit luminance
                l_limit = (selection(3)/100 + 1) .* obj.backgroundIntensity .* 255;
                u_limit = (selection(4)/100 + 1) .* obj.backgroundIntensity .* 255;
                rect_site = selection(6); % Site for rectifier
                
                % Pre-subunit rectification
                if rect_site == 0 && sum(selection(3:4) == Inf) == 0 % Inf = nonrectified condition
                    tmp(tmp < l_limit) = l_limit;
                    tmp(tmp > u_limit) = u_limit;
                    disp(strcat(mat2str([selection(3) selection(4)]),'% pre-subunit rectification applied...'))
                end
                
                % Apply subunit blur
                if subunitBlur_arcmin > 0
                    tmp = imgaussfilt(tmp,subunitBlur_arcmin);
                    disp(strcat(mat2str(selection(2)),'um subunit blur applied...'))
                end
                
                %%% After blur (before rectification), apply binning
                binningSize = selection(7);
                if binningSize > 0
                    
                    % Relative to background
                    background = obj.backgroundIntensity .* 255;
                    binningSize_8bit = background .* binningSize/100;
                    
                    % Symmetrically scale from 0% contrast
                    b_up = background:binningSize_8bit:255*2;
                    b_bdown = background:-binningSize_8bit:-255;
                    bins = [fliplr(b_bdown) b_up]; 
                    
                    % Remove bins (beyond 8-bit intensities)
                    bins(bins < 0) = [];
                    bins(bins > 255) = [];
                    bins = unique(bins);
                    
                    % Round to bins
                    tmp = interp1(bins,bins,tmp,'nearest');
                    disp(strcat(mat2str(selection(7)),'% contrast binning (posterization) applied...'))
                end
                
                % Post-subunit rectification
                if rect_site == 1 && sum(selection(3:4) == Inf) == 0 % Rectification = inf = nonrectified condition
                    tmp(tmp < l_limit) = l_limit;
                    tmp(tmp > u_limit) = u_limit;
                    disp(strcat(mat2str([selection(3) selection(4)]),'% post-subunit rectification applied...'))
                end
                
                % Apply RGC blur
                if rgcBlur_arcmin > 0
                    tmp = imgaussfilt(tmp,rgcBlur_arcmin);
                    disp(strcat(mat2str(selection(3)),'um subunit blur applied...'))
                end
                
                % (If posterization is used): how may unique luminance values are left?
                if binningSize > 0
                    U = unique(tmp);
                    U = round((U ./ background - 1) .* 100); % Convert to percent contrast
                    U(isnan(U)) = [];
                    disp(strcat(mat2str(length(U)), ' unique luminances present:',mat2str(U)))
                end
            end
            disp('...')
            
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
                aperture.size = [max(canvasSize) max(canvasSize)].*2;
                mask = stage.core.Mask.createCircularAperture(maskDiameterPix/max(canvasSize)/2, 1024); %circular aperture
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
