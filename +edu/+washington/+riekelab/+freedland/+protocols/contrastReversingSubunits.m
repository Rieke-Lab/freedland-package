% Searches for subunits within a center RF using contrast-reversing stimuli
classdef contrastReversingSubunits < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime     = 250   % ms
        stimTime    = 2000  % ms
        tailTime    = 250   % ms

        % Disk sizing and properties
        centerDiameter      = 200   % in um
        subunitDiameter     = 40    % in um
        primaryContrast     = 0.2   % contrast presented to center RF (0 to 1)
        dotContrast         = 0.9   % contrast presented to subunit-searching dot (0 to 1)
        temporalFrequency   = 4     % Hz
        
        % Additional options
        samplingOverlap     = 0.5;  % fractional overlap for subunit-searching dot (0 - 1)
        randomize           = false;    % randomize order of epochs
        backgroundIntensity = 0.168     % light intensity of background (0-1)
        onlineAnalysis      = 'extracellular'
        numberOfAverages    = uint16(1) % number of repeats to queue
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        centerDisk
        opposingDisk
        cartesianCoordinates
        polarCoordinates
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
            
            % Graphical figure
            obj.showFigure('edu.washington.riekelab.freedland.figures.contrastReversingSubunitsFigure',...
                obj.rig.getDevice(obj.amp),'temporalFrequency',obj.temporalFrequency,...
                'preTime',obj.preTime,'stimTime',obj.stimTime,...
                'monitorSampleRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'),...
                'centerRadius',obj.centerDiameter/2,'subunitRadius',obj.subunitDiameter/2); % in um

            % Convert units from microns to pixels
            canvasSize = fliplr(obj.rig.getDevice('Stage').getCanvasSize());
            centerRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.centerDiameter/2,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            subunitRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.subunitDiameter/2,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            
            % Create pixel space
            [xx,yy] = meshgrid(1:canvasSize(2),1:canvasSize(1));
            r = sqrt((xx - canvasSize(2)/2).^2 + (yy - canvasSize(1)/2).^2); 
            obj.centerDisk = r < centerRadiusPix; % Center RF mask
       
            % Tile subunit-searching dot locations
            radialDistance = 2 * subunitRadiusPix * (1-obj.samplingOverlap); % Euclidean distance
            radialCoordinates = 0:radialDistance:centerRadiusPix; % All sampling distances from center
            obj.opposingDisk = cell(1000,1);
            obj.cartesianCoordinates = zeros(1000,2);
            obj.polarCoordinates = zeros(1000,2);
            tempCounter = 1;
            for a = 1:length(radialCoordinates)
                side = radialCoordinates(a)/radialDistance; % Relative distance from center (in units of radii)
                angleSweep = acos(((side)^2 + (side).^2 - 1) / (2*side^2)); % Law of cosines
                angleCoordinates = 0:angleSweep:2*pi; % All possible angles
                for b = 1:length(angleCoordinates)
                    % Identify coordinates
                    obj.polarCoordinates(tempCounter,:) = [edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                        radialCoordinates(a),obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'pix2um'),...
                        rad2deg(angleCoordinates(b))];
                    
                    obj.cartesianCoordinates(tempCounter,1) = canvasSize(2)/2 + radialCoordinates(a) .* cos(angleCoordinates(b)); % x
                    obj.cartesianCoordinates(tempCounter,2) = canvasSize(1)/2 - radialCoordinates(a) .* sin(angleCoordinates(b)); % y
                    
                    % Make mask
                    r_subunit = sqrt((xx - obj.cartesianCoordinates(tempCounter,1)).^2 + (yy - obj.cartesianCoordinates(tempCounter,2)).^2) <= subunitRadiusPix; 
                    obj.opposingDisk{tempCounter,1} = r_subunit; % Subunit-searching dot
                    tempCounter = tempCounter + 1;
                end
            end

            % Remove zeros
            obj.opposingDisk(sum(obj.cartesianCoordinates,2) == 0,:) = [];
            obj.polarCoordinates(sum(obj.cartesianCoordinates,2) == 0,:) = [];
            obj.cartesianCoordinates(sum(obj.cartesianCoordinates,2) == 0,:) = [];

            % Check uniqueness
            [obj.cartesianCoordinates,adj] = unique(obj.cartesianCoordinates,'rows','stable');
            obj.polarCoordinates = obj.polarCoordinates(adj,:);
            obj.opposingDisk = obj.opposingDisk(adj,:);
            
            % Adjust cartesian coordinates for export
            obj.cartesianCoordinates = obj.cartesianCoordinates - repmat(fliplr(canvasSize/2),size(adj,1),1); % Center about zero
            obj.cartesianCoordinates = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.cartesianCoordinates,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'pix2um'); % Convert to microns
            
            % Estimate stimulus time
            disp(strcat('Unique locations to probe: ',mat2str(size(obj.opposingDisk,1))));
            expectedTime = (obj.preTime + obj.stimTime + obj.tailTime + 1.5*1000) / (1000*60); % +1.5 sec rig delay, in min
            disp(strcat('Approx. stimulus time: ',mat2str(round(size(obj.opposingDisk,1)*expectedTime*obj.numberOfAverages)),'min'));

            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(size(obj.cartesianCoordinates,1));
            else
                obj.order = 1:size(obj.cartesianCoordinates,1);
            end
        end
        
        function prepareEpoch(obj, epoch)

            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            % Add coordinates
            epoch.addParameter('subunitCoordinates_cartesianMicrons',obj.cartesianCoordinates(obj.order(obj.counter+1),:)); % About center (for figure)
            epoch.addParameter('subunitCoordinates_polarMicrons',obj.polarCoordinates(obj.order(obj.counter+1),:)); % For further analysis
            
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
            specificDisks = cat(3,obj.centerDisk,obj.opposingDisk{obj.order(obj.counter+1),1});

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
                    spContrast  = obj.primaryContrast; % Multiplier on square wave
                else
                    grate.phase = 0;
                    spContrast  = obj.dotContrast; % Multiplier on square wave
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