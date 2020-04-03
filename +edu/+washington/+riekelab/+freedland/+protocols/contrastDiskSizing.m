% Uses cell's RF to find NLs between center/surround
classdef contrastDiskSizing < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms

        % Disk sizing and properties
        diskRadii = [50 100]; % in um. automatically places a disk around the outer edge of the monitor
        contrast = 0.9 % relative to mean (0-1)
        temporalFrequency = 4 % Hz
        
        % Additional options
        splitField = false 
        rotation = 0;  % deg (0 to 180)
        backgroundIntensity = 0.168 % (0-1)
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(1) % number of epochs to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   
        diskRadii_PIX
        monitorFrameRate
        micronsPerPixel
        monitorSize
        analysisFigure
        disks
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
            
            % M. Turner's online analysis for F1/F2 frequencies
            if ~strcmp(obj.onlineAnalysis,'none')
                % custom figure handler
                if isempty(obj.analysisFigure) || ~isvalid(obj.analysisFigure)
                    obj.analysisFigure = obj.showFigure('symphonyui.builtin.figures.CustomFigure', @obj.F1F2_PSTH);
                    f = obj.analysisFigure.getFigureHandle();
                    set(f, 'Name', 'Cycle avg PSTH');
                    obj.analysisFigure.userData.runningTrace = 0;
                    obj.analysisFigure.userData.axesHandle = axes('Parent', f);
                else
                    obj.analysisFigure.userData.runningTrace = 0;
                end
            end
            
            % Pull variables
            obj.micronsPerPixel = obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.monitorSize = obj.rig.getDevice('Stage').getCanvasSize();
            obj.monitorSize = fliplr(obj.monitorSize); % Adjust to [height, width]
            obj.monitorFrameRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');

            % Identify disk radii
            obj.diskRadii_PIX = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.diskRadii,obj.micronsPerPixel,'UM2PIX'));
            obj.diskRadii_PIX = [0, obj.diskRadii_PIX, round(max(obj.monitorSize)/2)];
            
            % Make disks
            [xx,yy] = meshgrid(1:obj.monitorSize(2),1:obj.monitorSize(1));
            r = sqrt((xx - obj.monitorSize(2)/2).^2 + (yy - obj.monitorSize(1)/2).^2); 
            
            obj.disks = zeros(obj.monitorSize(1),obj.monitorSize(2),length(obj.diskRadii_PIX)-1);
            for a = 1:length(obj.diskRadii_PIX)-1
                obj.disks(:,:,a) = r > obj.diskRadii_PIX(a) & r < obj.diskRadii_PIX(a+1);
            end
            
            if obj.splitField == true;
                
                % Define theta space
                th = atan((xx - obj.monitorSize(2)/2) ./ (yy - obj.monitorSize(1)/2));
                th = abs(th-pi/2);              
                nonsmooth = find(diff(th) > pi/2,1);
                th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
                th = rad2deg(th);
                
                % Region for split field
                actualRotation = 180 - obj.rotation;
                overmask = th > actualRotation & th < actualRotation+180;
                invOvermask = overmask == 0;
                
                % Split masks
                set1 = obj.disks .* repmat(overmask,1,1,size(obj.disks,3));
                set2 = obj.disks .* repmat(invOvermask,1,1,size(obj.disks,3));
                obj.disks = cat(3,set1,set2);
                
                % Ensure adjacent disks are contrast reversing
                diskTotal = size(set1,3);
                partialOrder = [1 1+diskTotal];
                order = []; % Not large enough to bother pre-allocating.
                counter = 1;
                for a = 1:diskTotal
                    specificOrder = partialOrder + (counter-1)*2;
                    if mod(a,2) == 1
                        order = [order specificOrder];
                    else
                        order = [order  fliplr(specificOrder+1)];
                        counter = counter + 1;
                    end
                end

                % Order disks by contrast
                obj.disks = obj.disks(:,:,order);
            end
        end
        
        function F1F2_PSTH(obj, ~, epoch) %online analysis function

            % Function by M. Turner
            response = epoch.getResponse(obj.rig.getDevice(obj.amp));
            quantities = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            
            axesHandle = obj.analysisFigure.userData.axesHandle;
            runningTrace = obj.analysisFigure.userData.runningTrace;
            
            if strcmp(obj.onlineAnalysis,'extracellular') %spike recording
                filterSigma = (20/1000)*sampleRate; %msec -> dataPts
                newFilt = normpdf(1:10*filterSigma,10*filterSigma/2,filterSigma);
                res = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(quantities,[],sampleRate);
                epochResponseTrace = zeros(size(quantities));
                epochResponseTrace(res.sp) = 1; %spike binary
                epochResponseTrace = sampleRate*conv(epochResponseTrace,newFilt,'same'); %inst firing rate
            else %intracellular - Vclamp
                epochResponseTrace = quantities-mean(quantities(1:sampleRate*obj.preTime/1000)); %baseline
                if strcmp(obj.onlineAnalysis,'exc') %measuring exc
                    epochResponseTrace = epochResponseTrace./(-60-0); %conductance (nS), ballpark
                elseif strcmp(obj.onlineAnalysis,'inh') %measuring inh
                    epochResponseTrace = epochResponseTrace./(0-(-60)); %conductance (nS), ballpark
                end
            end
            
            noCycles = floor(obj.temporalFrequency*obj.stimTime/1000);
            period = (1/obj.temporalFrequency)*sampleRate; %data points
            epochResponseTrace(1:(sampleRate*obj.preTime/1000)) = []; %cut out prePts
            cycleAvgResp = 0;
            for c = 1:noCycles
                cycleAvgResp = cycleAvgResp + epochResponseTrace((c-1)*period+1:c*period);
            end
            cycleAvgResp = cycleAvgResp./noCycles;
            timeVector = (1:length(cycleAvgResp))./sampleRate; %sec
            runningTrace = runningTrace + cycleAvgResp;
            cla(axesHandle);
            h = line(timeVector, runningTrace./obj.numEpochsCompleted, 'Parent', axesHandle);
            set(h,'Color',[0 0 0],'LineWidth',2);
            xlabel(axesHandle,'Time (s)')
            title(axesHandle,'Running cycle average...')
            if strcmp(obj.onlineAnalysis,'extracellular')
                ylabel(axesHandle,'Spike rate (Hz)')
            else
                ylabel(axesHandle,'Resp (nS)')
            end
            obj.analysisFigure.userData.runningTrace = runningTrace;
        end
        
        function p = createPresentation(obj)
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity

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