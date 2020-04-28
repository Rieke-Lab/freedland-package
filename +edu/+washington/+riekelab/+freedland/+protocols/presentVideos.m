% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef presentVideos < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % RF field information
        rfSigmaCenter = 30; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.

        randomize = true; % whether to randomize movies shown
        rawMovieFrequency = 3; % How often should we show the original movie amongst experimental movies? (i.e, 3 = every 3 experimental movies, we show 1 raw movie).

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
            
            % General directory
            obj.directory = 'Documents/freedland-package/+edu/+washington/+riekelab/+freedland/+movies';
            D = dir(obj.directory);
            
            % Find correct folder
            seedName = [mat2str(obj.rfSigmaCenter),'_',mat2str(obj.rfSigmaSurround)];
            for a = 1:length(D)
                if strfind(D(a).name,seedName) > 0
                    folderName = D(a).name;
                end
            end
            
            % Correct directory, find movies
            obj.directory = [obj.directory,'/',folderName];
            D = dir(obj.directory);
            
            rawMovies = cell(size(D,1),1);
            testMovies = cell(size(D,1),1);
            for a = 1:length(D)
                if sum(strfind(D(a).name,'raw')) > 0
                    rawMovies{a,1} = D(a).name;
                elseif sum(strfind(D(a).name,'.mp4')) > 0
                    testMovies{a,1} = D(a).name;
                end
            end
            testMovies = testMovies(~cellfun(@isempty, testMovies(:,1)), :);
            rawMovies = rawMovies(~cellfun(@isempty, rawMovies(:,1)), :);
            
            % Seed in original files
            obj.movieFilenames = {};
            frequencyCheck = 0;
            rawMovieCounter = 1;
            for a = 1:size(testMovies,1)
                if frequencyCheck == 0
                    obj.movieFilenames = [obj.movieFilenames; {rawMovies{rawMovieCounter+1,1}}];
                    rawMovieCounter = mod(rawMovieCounter + 1,size(rawMovies,1));
                end
                
                obj.movieFilenames = [obj.movieFilenames; {testMovies{a,1}}];
                frequencyCheck = mod(frequencyCheck + 1,obj.rawMovieFrequency);
            end
            
            % Find background intensity
            [~, ~, ~, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(81, 1,...
                    'amplification', 1,'mirroring', false);
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            
            obj.sequence = repmat(1:size(obj.movieFilenames,1),1,obj.numberOfAverages);
            
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
            obj.counter = 1;

            % For identifying a good empirical rawMovieFrequency
            movieTimings(obj);
        
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = 6;
            
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
            p = stage.core.Presentation(6);

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
        
        function movieTimings(obj)
            videoLength = 6; % seconds per movie
            A = strfind(obj.movieFilenames,'raw');
            A = find(~cellfun(@isempty,A));
            
            % How often raw movie is played (for MATLAB 2016)
            B = [];
            for a = 1:length(A)
                B = [B,find(obj.sequence == A(a))];
            end
            
            C = diff(sort(B) .* videoLength / 60);
            D = strcat('Average time between raw movies:',mat2str(round(mean(C),2)),'minutes.');
            disp(D)
            
            D = strcat('Approx. total stimulus time (+33% rig delay):',mat2str(length(obj.sequence)* videoLength /60 * 1.33),'minutes.');
            disp(D)
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.sequence);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.sequence);
        end
    end
end