% Contrast-reversing spots that produce varying temporal shifts throughout
% the receptive field center. Eclectic set of 12 stimuli.
classdef julianSpatialBowlOfMojo < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        % Basic parameters
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms

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
        centerDisk
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
            
            obj.centerDisk = (r <= centerRadius_px);
            
            % Build different gratings conditions
            gratingWedges = [2; 8];
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
            
            % Define stimuli in Julian's spatial bowl of mojo
            obj.combinatoricDivisions = [...
                0	0	0	0	0
                1	0	0	0	1
                2	0	0	0	2
                0	2	2	2	0
                0	0	1	1	1
                0	0	2	2	2
                2	2	0	2	0
                1	1	0	1	0
                0	0	1	0	0
                0	0	2	0	0
                1	0	1	0	1
                2	3	2	3	2];
            
            % Interweave both grating conditions to build "distorted" stimuli
            obj.masks = cell(size(obj.combinatoricDivisions,1),1);
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

            % Add contrast reversing gratings
            for a = 1:2
                grate           = stage.builtin.stimuli.Grating('square'); % Square wave grating
                grate.size      = canvasSize;
                grate.position  = canvasSize / 2;
                grate.spatialFreq = 1e-5;
                grate.color     = 2 * obj.backgroundIntensity; % Amplitude of square wave
                grate.contrast  = obj.contrast; % Multiplier on square wave
                
                % Define in-phase and out-of-phase mask
                specificMask = obj.masks{obj.sequence(obj.counter+1),1};
                if a == 2
                    specificMask = mod(specificMask + 1,2) .* obj.centerDisk;
                end
                grateShape      = uint8(specificMask*255);
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