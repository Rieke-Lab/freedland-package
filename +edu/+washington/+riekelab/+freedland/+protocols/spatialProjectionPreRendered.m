% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef spatialProjectionPreRendered < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        rfSigmaCenter       = 70; % in um
        rfSigmaSurround     = 170; % in um
        diskRadii           = [100,200,300]; % in um
        naturalImages       = [5,81];
        
        rawTrajectoryFrequency = 5; % how often to see original image
        randomize = true

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        backgroundIntensity
        sequence
        counter
        imageMatrix
        movieFilenames
        directory
        totalRuns
        videoLength = 6; % in sec
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
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',1);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            fileSettings = strcat(mat2str(obj.rfSigmaCenter),'_',mat2str(obj.rfSigmaSurround),'_',mat2str(obj.diskRadii),'\');
                        
            obj.directory = strcat('Documents\freedland-package\+edu\+washington\+riekelab\+freedland\+movies\',fileSettings);
            D = dir(obj.directory);
            
            % Pull relevant movie filenames
            filenames = [];
            rawFilenames = [];
            for b = 1:length(obj.naturalImages)
                imgName = strcat('img',mat2str(obj.naturalImages(b)));
                for a = 1:size(D)
                    if sum(strfind(D(a).name,imgName)) > 0
                        if sum(strfind(D(a).name,'raw')) == 0
                            filenames = [filenames; {D(a).name}];
                        else
                            rawFilenames = [rawFilenames; {D(a).name}];
                        end
                    end
                end
            end
            
            % Seed in original files
            obj.movieFilenames = [];
            frequencyCheck = 0;
            rawMovieCounter = 1;
            for a = 1:size(filenames,1)
                if frequencyCheck == 0
                    obj.movieFilenames = [obj.movieFilenames; {rawFilenames{rawMovieCounter,1}}];
                    rawMovieCounter = mod(rawMovieCounter + 1,size(rawFilenames,1)) + 1;
                end
                
                obj.movieFilenames = [obj.movieFilenames; {filenames{a,1}}];
                frequencyCheck = mod(frequencyCheck + 1,obj.rawTrajectoryFrequency);
            end
            

            % Find background intensity
            [~, ~, ~, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.naturalImages(1), 1,...
                    'amplification', 1,'mirroring', false);
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            
            obj.sequence = repmat(1:size(obj.movieFilenames,1),1,obj.numberOfAverages);
            
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
            obj.counter = 1;
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = obj.videoLength;
            
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
            p = stage.core.Presentation(obj.videoLength);

            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            selection = obj.sequence(obj.counter);
            f = obj.movieFilenames{selection,1}; % filename for relevant movie

            % Prep to display image
            scene = stage.builtin.stimuli.Movie(fullfile(obj.directory,f));
            scene.size = [canvasSize(1),canvasSize(2)];
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            p.addStimulus(scene);
            
            obj.counter = obj.counter + 1;
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.sequence);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.sequence);
        end
    end
end