% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFMetamer < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime = 250 % in ms
        stimTime = 5500 % in ms
        tailTime = 250 % in ms
        
        % Natural image trajectory
        referenceImage = 81; % image to build frankenstein metamer for.
        
        % RF field information
        rfSigmaCenter = 30; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.
        
        % Disk placement
        disks = 3; % number of disks to replace image with.
        overrideRadii = [0 0.75 2 3]; % only takes effect if any value is >0. Arranges disks in any distribution depending on coordinate system. For info on RF coordinate system, see RFConversion function in code.
        overrideCoordinate = 'RF'; % type of coordinates to measure disk radii.
        xSliceFrequency = 1; % how many radial slices to cut between 0 and 90 degrees.
        ySliceFrequency = 1; % how many radial slices to cut between 90 and 180 degrees.
        disksIgnoreCut = [0 3]; % starting from the center disk and moving outwards, how many disks should we NOT cut (keep circular)?

        % Disk type
        meanDisks = [1 2 3]; % starting from the center disk and moving outwards, which disks should be averaged?
        naturalDisks = [0 0];  % starting from the center disk and moving outwards, which disks should remain a natural image?
        backgroundDisks = [0 0]; % starting from the center disk and moving outwards, which disks should be left at background intensity?
        switchDisks = [0 0]; % starting from the center disk and moving outwards, which disks should switch intensity based on fixation?
        meanIntegration = 'gaussian'; % type of linear integration
        
        % Additional parameters
        numberOfDistinctMovies = 5; % number of distinct movies to compare
        randomize = true;
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(8) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        meanIntegrationType = symphonyui.core.PropertyType('char', 'row', {'uniform','gaussian'})
        backgroundIntensity
        sequence
        counter
        imageMatrix
        movieFilenames
        directory
    end

    methods
        
        function didSetRig(obj)
            didSetRig@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            [obj.amp, obj.ampType] = obj.createDeviceNamesProperty('Amp');
        end

        function prepareRun(obj)

            prepareRun@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj);
            
            % Must have prerender set to avoid lagging through the replacement trajectory.     
            if obj.rig.getDevice('Stage').getConfigurationSetting('prerender') == 0
                error('Must have prerender set') 
            end

            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',2);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
                        
            D = dir('+edu/+washington/+riekelab/+freedland/+movies');
            obj.movieFilenames = [];
            replacementMovies = [];
            for a = 1:size(D,1)
                A = D(a).name;
                
                % Only select relevant videos
                if contains(A,string(obj.referenceImage)) && ...
                        contains(A,string(obj.rfSigmaCenter)) && ...
                        contains(A,string(obj.rfSigmaSurround))
                    if contains(A,'ref') % Reference trajectory, must include
                        obj.movieFilenames = [obj.movieFilenames;{A}];
                    elseif contains(A,'rep') % Replacement trajectories
                        replacementMovies = [replacementMovies;{A}];
                    end
                end
            end
            replacementMovies = replacementMovies(1:obj.numberOfDistinctMovies,:);
            obj.movieFilenames = [obj.movieFilenames;replacementMovies];

            obj.directory = '+edu/+washington/+riekelab/+freedland/+movies/';
            uiopen(fullfile(obj.directory,'ref81.mp4'))
            
            raw = double(ref81);
            obj.backgroundIntensity = raw(1,1,1,1) / 255; % set manually in videos in top left corner
            
            obj.sequence = repelem(1:obj.numberOfDistinctMovies,obj.numberOfAverages);
            
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
            obj.counter = 1;
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('movieName',obj.movieFilenames{obj.sequence(obj.counter),1});

            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);

            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            selection = obj.sequence(obj.counter);
            f = obj.movieFilenames{selection,1}; % filename for relevant movie

            % Prep to display image
            scene = stage.builtin.stimuli.Movie(fullfile(obj.directory,f));
            scene.size = [size(obj.imageMatrix,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(obj.imageMatrix,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            p.addStimulus(scene);
            
            obj.counter = obj.counter + 1;
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages*obj.numberOfDistinctMovies;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages*obj.numberOfDistinctMovies;
        end
    end
end