classdef contrastReversingGratingsCenterSurround < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms
        centerBarWidth = [5 10 20 40 60 80] % um
        centerContrast = [0.05 0.1 0.2 0.3 0.6] % relative to mean (0-1)
        surroundBarWidth = [40 80 120 160 200 240 300] % um
        surroundContrast = [0.05 0.1 0.2 0.3 0.6] % relative to mean (0-1)
        surroundPhase = 'in-phase';
        temporalFrequency = 4 % Hz
        apertureDiameter = 300; % um
        rotation = 0; % deg
        backgroundIntensity = 0.168 % (0-1)
        randomizeOrder = true;
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        surroundPhaseType = symphonyui.core.PropertyType('char', 'row', {'in-phase', 'out-of-phase'})
        sequence
        counter
        xAxis
        li
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
            
            % Arrange experiment conditions
            counter2 = 1;
            obj.sequence = zeros(length(obj.centerContrast).*length(obj.centerBarWidth).*length(obj.surroundContrast).*length(obj.surroundBarWidth),4);
            for a = 1:length(obj.centerContrast)
                for b = 1:length(obj.centerBarWidth)
                    for c = 1:length(obj.surroundContrast)
                        for d = 1:length(obj.surroundBarWidth)
                            obj.sequence(counter2,:) = [obj.centerContrast(a), obj.centerBarWidth(b), obj.surroundContrast(c), obj.surroundBarWidth(d)];
                            counter2 = counter2 + 1;
                        end
                    end
                end
            end
            
            % Identify experiment condition with most settings
            l = [];
            for a = 1:size(obj.sequence,2)
                l(a) = length(unique(obj.sequence(:,a)));
            end
            [~,obj.li] = max(l);
            obj.xAxis = unique(obj.sequence(:,obj.li)); % Define x-axis for figure.
            
            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',1);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            %%% M. Turner's analysis figure
            if length(obj.xAxis) > 1
                colors = edu.washington.riekelab.turner.utils.pmkmp(length(obj.xAxis),'CubicYF');
            else
                colors = [0 0 0];
            end
            if ~strcmp(obj.onlineAnalysis,'none')
                % Custom figure handler
                if isempty(obj.analysisFigure) || ~isvalid(obj.analysisFigure)
                    obj.analysisFigure = obj.showFigure('symphonyui.builtin.figures.CustomFigure', @obj.CRGanalysis);
                    f = obj.analysisFigure.getFigureHandle();
                    set(f, 'Name', 'CRGs');
                    obj.analysisFigure.userData.trialCounts = zeros(size(obj.xAxis));
                    obj.analysisFigure.userData.F1 = zeros(size(obj.xAxis));
                    obj.analysisFigure.userData.F2 = zeros(size(obj.xAxis));
                    obj.analysisFigure.userData.axesHandle = axes('Parent', f);
                else
                    obj.analysisFigure.userData.trialCounts = zeros(size(obj.xAxis));
                    obj.analysisFigure.userData.F1 = zeros(size(obj.xAxis));
                    obj.analysisFigure.userData.F2 = zeros(size(obj.xAxis));
                end
            end
            %%%
            
            % Arrange experiment conditions
            counter2 = 1;
            obj.sequence = zeros(length(obj.centerContrast).*length(obj.centerBarWidth).*length(obj.surroundContrast).*length(obj.surroundBarWidth),4);
            for a = 1:length(obj.centerContrast)
                for b = 1:length(obj.centerBarWidth)
                    for c = 1:length(obj.surroundContrast)
                        for d = 1:length(obj.surroundBarWidth)
                            obj.sequence(counter2,:) = [obj.centerContrast(a), obj.centerBarWidth(b), obj.surroundContrast(c), obj.surroundBarWidth(d)];
                            counter2 = counter2 + 1;
                        end
                    end
                end
            end
   
            % Randomize order
            if obj.randomizeOrder == true
                obj.sequence = obj.sequence(randperm(size(obj.sequence,1)),:);
            end

            obj.counter = 1;
            obj.sequence = repmat(obj.sequence,obj.numberOfAverages,1);
            
            %%% Quick bug fix: very last epoch won't play, so we double last
            %%% condition and tell Symphony to present one less epoch
            obj.sequence = [obj.sequence; obj.sequence(end,:)];
        end
        
        function CRGanalysis(obj, ~, epoch) % Online analysis function by M. Turner
            response = epoch.getResponse(obj.rig.getDevice(obj.amp));
            epochResponseTrace = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            
            axesHandle = obj.analysisFigure.userData.axesHandle;
            trialCounts = obj.analysisFigure.userData.trialCounts;
            F1 = obj.analysisFigure.userData.F1;
            F2 = obj.analysisFigure.userData.F2;
            
            if strcmp(obj.onlineAnalysis,'extracellular') %spike recording
                %take (prePts+1:prePts+stimPts)
                epochResponseTrace = epochResponseTrace((sampleRate*obj.preTime/1000)+1:(sampleRate*(obj.preTime + obj.stimTime)/1000));
                %count spikes
                S = edu.washington.riekelab.turner.utils.spikeDetectorOnline(epochResponseTrace);
                epochResponseTrace = zeros(size(epochResponseTrace));
                epochResponseTrace(S.sp) = 1; %spike binary
                
            else %intracellular - Vclamp
                epochResponseTrace = epochResponseTrace-mean(epochResponseTrace(1:sampleRate*obj.preTime/1000)); %baseline
                %take (prePts+1:prePts+stimPts)
                epochResponseTrace = epochResponseTrace((sampleRate*obj.preTime/1000)+1:(sampleRate*(obj.preTime + obj.stimTime)/1000));
            end

            L = length(epochResponseTrace); %length of signal, datapoints
            X = abs(fft(epochResponseTrace));
            X = X(1:L/2);
            f = sampleRate*(0:L/2-1)/L; %freq - hz
            [~, F1ind] = min(abs(f-obj.temporalFrequency)); %find index of F1 and F2 frequencies
            [~, F2ind] = min(abs(f-2*obj.temporalFrequency));

            F1power = 2*X(F1ind); %pA^2/Hz for current rec, (spikes/sec)^2/Hz for spike rate
            F2power = 2*X(F2ind); %double b/c of symmetry about zero
            
            barInd = find(obj.sequence(obj.counter-1,obj.li) == obj.xAxis); % -1 because we already ran update
            trialCounts(barInd) = trialCounts(barInd) + 1;
            F1(barInd) = F1(barInd) + F1power;
            F2(barInd) = F2(barInd) + F2power;
            
            cla(axesHandle);
            h1 = line(obj.xAxis, F1./trialCounts, 'Parent', axesHandle);
            set(h1,'Color','g','LineWidth',2,'Marker','o');
            h2 = line(obj.xAxis, F2./trialCounts, 'Parent', axesHandle);
            set(h2,'Color','r','LineWidth',2,'Marker','o');
            hl = legend(axesHandle,{'F1','F2'});
            xlabel(axesHandle,'Bar width (um)')
            ylabel(axesHandle,'Amplitude')

            obj.analysisFigure.userData.trialCounts = trialCounts;
            obj.analysisFigure.userData.F1 = F1;
            obj.analysisFigure.userData.F2 = F2;
        end
             
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('currentCenterContrast', obj.sequence(obj.counter,1));
            epoch.addParameter('currentCenterBarWidth', obj.sequence(obj.counter,2));
            epoch.addParameter('currentSurroundContrast', obj.sequence(obj.counter,3));
            epoch.addParameter('currentSurroundBarWidth', obj.sequence(obj.counter,4));
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