% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef spatialProjectionPreRendered < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        rfSigmaCenter       = 70;
        rfSigmaSurround     = 170;
        diskRadii           = [100,200,300]; % in um
        naturalImages       = [5,81];

        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        rfSigmaCenterType = symphonyui.core.PropertyType('char', 'row', {'30','40','50','70','100'}) 
        rfSigmaSurroundType = symphonyui.core.PropertyType('char', 'row', {'100','130','150','180','210'}) 
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
                        
            obj.directory = 'Documents\freedland-package\+edu\+washington\+riekelab\+freedland\+movies';
            D = dir(obj.directory);
            
            obj.movieFilenames = [];
            frequencyCheck = 0;
            
            for j = 1:length(obj.referenceImage) % Individual image
                
                % Find reference movie
                refMovie = strcat('img',mat2str(obj.referenceImage(j)),'_raw');
                for a = 1:size(D,1)
                    A = D(a).name;

                    if sum(strfind(A,refMovie)) > 0 % Normal trajectory
                        movieName = A;
                    end
                end
                
                % Find relevant diskArray videos
                for b = 1:length(obj.experimentsDiskArray)

                    specificFile = strcat('img',mat2str(obj.referenceImage(j)),'_',obj.rfSigmaCenter,...
                        '_',obj.rfSigmaSurround,'_',mat2str(obj.experimentsDiskArray(b)),'_projection');
                    specificFileBlur = strcat('img',mat2str(obj.referenceImage(j)),'_',obj.rfSigmaCenter,...
                        '_',obj.rfSigmaSurround,'_',mat2str(obj.experimentsDiskArray(b)),'_blur__projection');

                    replacementMovies = [];
                    
                    for a = 1:size(D,1)
                        A = D(a).name;

                        if sum(strfind(A,specificFile)) > 0 % Replacement trajectory
                            
                            if mod(frequencyCheck,obj.rawMovieFrequency) == 0
                                replacementMovies = [replacementMovies;{movieName}];
                            end
                                
                            replacementMovies = [replacementMovies;{A}];
                            frequencyCheck = frequencyCheck + 1;
                        end
                        
                        if sum(obj.experimentsDiskArray(b) == obj.addBlur) > 0  % Add video with blur into the mix
                            if sum(strfind(A,specificFileBlur)) > 0
                                if mod(frequencyCheck,obj.rawMovieFrequency) == 0
                                    replacementMovies = [replacementMovies;{movieName}];
                                end

                                replacementMovies = [replacementMovies;{A}];
                                frequencyCheck = frequencyCheck + 1;
                            end
                        end
                            
                    end

                    obj.movieFilenames = [obj.movieFilenames;replacementMovies];
                end

                % Repeat for metamers
                for b = 1:length(obj.experimentsMetamer)

                    specificFile = strcat('img',mat2str(obj.referenceImage(j)),'_',obj.rfSigmaCenter,...
                        '_',obj.rfSigmaSurround,'_',mat2str(obj.experimentsMetamer(b)),'_metamer');
                    specificFileBlur = strcat('img',mat2str(obj.referenceImage(j)),'_',obj.rfSigmaCenter,...
                        '_',obj.rfSigmaSurround,'_',mat2str(obj.experimentsMetamer(b)),'_blur__metamer');

                    replacementMovies = [];
                    individualCounter = 0;
                    individualBlurCounter = 0;

                    for a = 1:size(D,1)
                        A = D(a).name;

                        if sum(strfind(A,specificFile)) > 0 % Replacement trajectory
                            
                            if individualCounter < obj.numberOfMetamers
                                
                                if mod(frequencyCheck,obj.rawMovieFrequency) == 0
                                    replacementMovies = [replacementMovies;{movieName}];
                                end
                                
                                replacementMovies = [replacementMovies;{A}];
                                frequencyCheck = frequencyCheck + 1;
                                individualCounter = individualCounter + 1;
                            end
                        end
                        
                        if sum(obj.experimentsMetamer(b) == obj.addBlur) > 0  % Add video with blur into the mix
                            if sum(strfind(A,specificFileBlur)) > 0
                                
                                if individualBlurCounter < obj.numberOfMetamers
                                    if mod(frequencyCheck,obj.rawMovieFrequency) == 0
                                        replacementMovies = [replacementMovies;{movieName}];
                                    end
                                    
                                    replacementMovies = [replacementMovies;{A}];
                                    frequencyCheck = frequencyCheck + 1;
                                    individualBlurCounter = individualBlurCounter + 1;
                                end

                            end
                        end
                    end

                    obj.movieFilenames = [obj.movieFilenames;replacementMovies];
                end
            end

            obj.totalRuns = size(obj.movieFilenames,1);

            % Find background intensity
            [~, ~, ~, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.referenceImage(1), 1,...
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
            movieTimings(obj,movieName);
        
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.totalStimulusTime) / 1e3;
            
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
            p = stage.core.Presentation(obj.totalStimulusTime * 1e-3);

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
        
        function movieTimings(obj,movieName)
            videoLength = 6; % seconds
            A = find(strcmp(obj.movieFilenames,movieName));
            
            B = [];
            for a = 1:length(A)
                B = [B,find(obj.sequence == A(a))];
            end
            
            C = diff(sort(B) .* videoLength / 60);
            D = strcat('Average time between raw movies:',mat2str(round(mean(C),2)),'minutes.');
            disp(D)
            
            D = strcat('Approximate total stimulus time:',mat2str(length(obj.sequence)* videoLength /60),'minutes.');
            disp(D)
        end

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages*obj.totalRuns;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages*obj.totalRuns;
        end
    end
end