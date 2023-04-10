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
            
            % Define dimensional space
            c_dim = length(obj.circularDivisions);
            
            % Define all sets of possible locations for spatial discontinuities
            obj.combinatoricDivisions = zeros(1,c_dim);
            for a = 1:c_dim
                tmp = nchoosek(1:c_dim,a);
                for b = 1:size(tmp,1)
                    tmp2 = zeros(1,c_dim);
                    tmp2(tmp(b,:)) = 1;
                    obj.combinatoricDivisions = [obj.combinatoricDivisions; tmp2];
                end
            end
            
            % Include additional conditions
            if obj.controlGratings ~= obj.replacedGratings
                obj.combinatoricDivisions = [obj.combinatoricDivisions; ...     % [original + rotated gratings]
                                             obj.combinatoricDivisions .* 2;... % [original + shuffled gratings]
                                             obj.combinatoricDivisions + 2];    % [mixed shuffled gratings]
            end
            obj.combinatoricDivisions = unique(obj.combinatoricDivisions,'rows');
            
            % Remove complements (already sampled via reversing contrasts)
            for a = 1:size(obj.combinatoricDivisions)
                tmp = obj.combinatoricDivisions(a,:);
            
                % Check for inverse cases (control gratings)
                if sum(tmp == 1 | tmp == 0) == length(tmp)
                    tmp2 = mod(tmp + 1,2);
                    [~, i] = ismember(obj.combinatoricDivisions,tmp2,'rows');
                    obj.combinatoricDivisions(logical(i),:) = NaN;
                end
            
                % Check for inverse cases (replaced gratings)
                if sum(tmp == 2 | tmp == 3) == length(tmp)
                    tmp2 = mod(tmp + 1,2) + 2;
                    [~, i] = ismember(obj.combinatoricDivisions,tmp2,'rows');
                    obj.combinatoricDivisions(logical(i),:) = NaN;
                end
            end
            obj.combinatoricDivisions = obj.combinatoricDivisions(~isnan(obj.combinatoricDivisions(:,1)),:);
            
            % Interweave both grating conditions to build "distorted" stimuli
            obj.masks = cell(size(obj.combinatoricDivisions,1),2);
            for a = 1:size(obj.combinatoricDivisions,1) % Each distorted stimulus
            
                maskMap = obj.combinatoricDivisions(a,:);
                tmp = zeros(size(r));
            
                for b = 1:length(maskMap) % Each outlying region
            
                    % Build annulus
                    innerDiameter = (b-1) / length(maskMap);
                    outerDiameter = (b) / length(maskMap);
                    maskArea = (r > (centerRadius_px .* innerDiameter) & r <= (centerRadius_px .* outerDiameter));
            
                    if maskMap(b) == 0 % Original grating
                        tmp = tmp + builtMasks{1,1} .* maskArea;
                    elseif maskMap(b) == 1 % Rotated original grating
                        tmp = tmp + builtMasks{1,2} .* maskArea;
                    elseif maskMap(b) == 2 % Shuffled replaced grating
                        tmp = tmp + builtMasks{2,1} .* maskArea;
                    elseif maskMap(b) == 3 % Rotated replaced grating
                        tmp = tmp + builtMasks{2,2} .* maskArea;
                    end
                end
                obj.masks{a,1} = tmp;
                obj.masks{a,2} = mod(tmp+1,2) .* (r <= centerRadius_px);
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