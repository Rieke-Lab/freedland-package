% Places two sets of gratings, separated by spatial discontinuity
classdef combinatoricDiscontinuities < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        % Basic parameters
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms

        % Regular gratings
        controlGratings = 2; % Number of wedges in baseline stimulus (2 = split-field spot)
        circularDivisions = [20, 40, 60, 80, 100]; % Location of spatial discontinuities (% of RF center)
        replacedGratings = 8; % Number of wedges for replaced stimulus
        replacedRegion = 'inner' % inner or outer center

        % Stimulus parameters
        centerDiameter = 300; % in um
        contrast = 0.9; % 0 - 1
        temporalFrequency = 4 % Hz. Set to 0 for spatial flash
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
        replacedRegionType = symphonyui.core.PropertyType('char', 'row', {'inner', 'outer'})
        masks
        sequence
        counter
        combinatoricDivisions
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

            % Starting parameters
            centerRadius_px = obj.rig.getDevice('Stage').um2pix(obj.centerDiameter)/2;
            canvasSize      = obj.rig.getDevice('Stage').getCanvasSize();

            % Define polar space
            [xx,yy] = meshgrid(1:canvasSize(1),1:canvasSize(2));
            r = sqrt((xx - canvasSize(1)/2).^2 + (yy - canvasSize(2)/2).^2);
            th = atan((xx - canvasSize(1)/2) ./ (yy - canvasSize(2)/2));
            th = abs(th-pi/2);
            nonsmooth = find(diff(th) > pi/2,1);
            th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
            th = rad2deg(th);
            
            % Build different gratings conditions
            gratingWedges = [obj.controlGratings; obj.replacedGratings];
            builtMasks = cell(length(gratingWedges),2);
            for iter = 1:2
                % Pre-allocate space
                for a = 1:2
                    builtMasks{iter,a} = zeros(size(r));
                end
            
                % Define wedges for each condition
                A = gratingWedges(iter);
                A(A == 0) = 1;
                theta = 0:360/A:360;
            
                % Build each grating
                for a = 1:length(theta) - 1
                    m = (th >= theta(a) & th < theta(a+1)) .* (r <= centerRadius_px);
                    builtMasks{iter,mod(a,2)+1} = builtMasks{iter,mod(a,2)+1} + m;
                end
            end
            
            % Define all sets of possible locations for spatial discontinuities
            obj.combinatoricDivisions = zeros(1,length(obj.circularDivisions)-1);
            for a = 1:length(obj.circularDivisions)-1
                tmp = nchoosek(1:length(obj.circularDivisions)-1,a);
                for b = 1:size(tmp,1)
                    tmp2 = zeros(1,length(obj.circularDivisions)-1);
                    tmp2(tmp(b,:)) = 1;
                    obj.combinatoricDivisions = [obj.combinatoricDivisions; tmp2];
                end
            end
            
            % Include additional conditions
            obj.combinatoricDivisions = [obj.combinatoricDivisions ones(size(obj.combinatoricDivisions,1),1)];
            obj.combinatoricDivisions = [zeros(1,size(obj.combinatoricDivisions,2)); obj.combinatoricDivisions];
            
            % Interweave both grating conditions to build "distorted" stimuli
            maskOrdering = [2 1];
            obj.masks = cell(size(obj.combinatoricDivisions,1),2);
            for a = 1:size(obj.masks,1) % Each distorted stimulus
                
                % Identify combinatoric pattern
                spatialDivisions = obj.circularDivisions(obj.combinatoricDivisions(a,:) == 1);
                tmp_stimulus = zeros(size(r));
                tmp_tracker = zeros(size(r));
                if ~isempty(spatialDivisions)
            
                    %%% Place gratings from inner center outwards
                    % First dimension regulates which type of gratings to add
                    % mod(b,2)+1 regulates whether to add in-phase or out-of-phase region
                    for b = 1:length(spatialDivisions)
                        innerArea = (r <= (centerRadius_px .* spatialDivisions(b)/100));  
            
                        if strcmp(obj.replacedRegion,'inner') 
                            % Replace gratings for all but outer center
                            if b < length(spatialDivisions)
                                tmp1 = builtMasks{maskOrdering(1),mod(b,2)+1} .* innerArea - tmp_tracker;
                            else
                                tmp1 = builtMasks{maskOrdering(2),mod(b,2)+1} .* innerArea - tmp_tracker;
                            end
                        elseif strcmp(obj.replacedRegion,'outer')
                            % Replace gratings for all but inner center
                            if b > 1
                                tmp1 = builtMasks{maskOrdering(1),mod(b,2)+1} .* innerArea - tmp_tracker;
                            else
                                tmp1 = builtMasks{maskOrdering(2),mod(b,2)+1} .* innerArea - tmp_tracker;
                            end
                        end
                        tmp1(tmp1 < 0) = 0;
                        tmp_tracker = tmp_tracker + innerArea;
                        tmp_stimulus = tmp_stimulus + tmp1;
                    end
                    obj.masks{a,1} = tmp_stimulus;
                else
                    obj.masks{a,1} = builtMasks{maskOrdering(1),1} .* (r <= (centerRadius_px));
                end
            end

            % Randomize order
            obj.sequence = 1:size(obj.masks,1);
            if obj.randomizeOrder == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end

            obj.counter = 0;
            obj.sequence = repmat(obj.sequence,1,obj.numberOfAverages);
        end

        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);

            epoch.addParameter('spatialDiscontinuities', obj.combinatoricDivisions(obj.sequence(obj.counter+1),:));
        end

        function p = createPresentation(obj)
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();

            specificMasks = obj.masks(obj.sequence(obj.counter+1),:);
            
            % Add contrast reversing gratings
            for a = 1:2
                grate           = stage.builtin.stimuli.Grating('square'); % Square wave grating
                grate.size      = canvasSize;
                grate.position  = canvasSize / 2;
                grate.spatialFreq = 1e-5;
                grate.color     = 2 * obj.backgroundIntensity; % Amplitude of square wave
                grate.contrast  = obj.contrast; % Multiplier on square wave
                grateShape      = uint8(specificMasks{a}*255);
                grateMask       = stage.core.Mask(grateShape);
                grate.contrast  = obj.contrast;
                grate.setMask(grateMask);
              
                % Contrasting gratings between each set of masks
                if a == 1
                    grate.phase = 180;
                else
                    grate.phase = 0;
                end

                p.addStimulus(grate); % Add grating to the presentation
                
                % Control contrast
                if obj.temporalFrequency > 0
                    grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                        @(state)getGrateContrast(obj, state.time - obj.preTime/1e3));
                    p.addController(grateContrast); % Add the controller
                end
                
                % Hide during pre & post
                grateVisible = stage.builtin.controllers.PropertyController(grate, 'visible', ...
                    @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
                p.addController(grateVisible);
            end
            obj.counter = obj.counter + 1;

            function c = getGrateContrast(obj, time)
                c = obj.contrast.*sin(2 * pi * obj.temporalFrequency * time);
            end
        end
 
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.sequence);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.sequence);
        end
        
    end
    
end