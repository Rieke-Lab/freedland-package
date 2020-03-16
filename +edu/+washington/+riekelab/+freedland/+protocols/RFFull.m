% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFFull < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Basics
        totalStimulusTime = 6000; % in ms
        referenceImage = 81; % image to build metamer for.
        
        % Arrangement of disks
        experimentsDiskArray = 1:15  % Number of experimental movies to show
        experimentsMetamer = [8,9,12,13,15]; % Pre-determined model number to show replacement over.
        
        % RF field information
        rfSigmaCenter = 30; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.

        addBlur = true; % apply gaussian blur over sharp edges.
        randomize = true; % whether to randomize movies shown
        numberOfMetamers = 3; % number of distinct metamers to show
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(8) % number of epochs to queue
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
                        
            obj.directory = 'Documents\freedland-package\+edu\+washington\+riekelab\+freedland\+movies';
            D = dir(obj.directory);
            
            obj.movieFilenames = [];
            % Find relevant diskArray videos
            for b = 1:length(obj.experimentsDiskArray)
                
                % Identify relevant filenames
                if obj.addBlur == false
                    specificFile = strcat('img',mat2str(obj.referenceImage),'_',mat2str(obj.rfSigmaCenter),...
                        '_',mat2str(obj.rfSigmaSurround),'_',mat2str(obj.experimentsDiskArray(b)),'_projection');
                else
                    specificFile = strcat('img',mat2str(obj.referenceImage),'_',mat2str(obj.rfSigmaCenter),...
                        '_',mat2str(obj.rfSigmaSurround),'_',mat2str(obj.experimentsDiskArray(b)),'_blur__projection');
                end
                
                refMovie = strcat('img',mat2str(obj.referenceImage),'_raw');
                replacementMovies = [];
                
                for a = 1:size(D,1)
                    A = D(a).name;

                    % Only select relevant videos
                    if sum(strfind(A,refMovie)) > 0 % Normal trajectory
                        obj.movieFilenames = [obj.movieFilenames;{A}];
                    end

                    if sum(strfind(A,specificFile)) > 0 % Replacement trajectory
                        replacementMovies = [replacementMovies;{A}];
                    end
                end

                obj.movieFilenames = [obj.movieFilenames;replacementMovies];
            end
            
            % Repeat for metamers
            for b = 1:length(obj.experimentsMetamer)
                
                % Identify relevant filenames
                if obj.addBlur == false
                    specificFile = strcat('img',mat2str(obj.referenceImage),'_',mat2str(obj.rfSigmaCenter),...
                        '_',mat2str(obj.rfSigmaSurround),'_',mat2str(obj.experimentsMetamer(b)),'_metamer');
                else
                    specificFile = strcat('img',mat2str(obj.referenceImage),'_',mat2str(obj.rfSigmaCenter),...
                        '_',mat2str(obj.rfSigmaSurround),'_',mat2str(obj.experimentsMetamer(b)),'_blur__metamer');
                end
                
                refMovie = strcat('img',mat2str(obj.referenceImage),'_raw');
                replacementMovies = [];
                
                for a = 1:size(D,1)
                    A = D(a).name;

                    % Only select relevant videos
                    if sum(strfind(A,refMovie)) > 0 % Normal trajectory
                        obj.movieFilenames = [obj.movieFilenames;{A}];
                    end

                    if sum(strfind(A,specificFile)) > 0 % Replacement trajectory
                        replacementMovies = [replacementMovies;{A}];
                    end
                end
                
                if length(replacementMovies) > obj.numberOfMetamers
                    replacementMovies = replacementMovies(1:obj.numberOfMetamers,:);
                elseif length(replacementMovies) < obj.numberOfMetamers
                    error('not enough metamers pre-generated.')
                end

                obj.movieFilenames = [obj.movieFilenames;replacementMovies];
            end

            obj.totalRuns = size(obj.movieFilenames,1);
                        
            % Find background intensity
            [~, ~, ~, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(obj.referenceImage, 1,...
                    'amplification', 1,'mirroring', false);
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            
            obj.sequence = repelem(1:size(obj.movieFilenames,1),obj.numberOfAverages); % +1 includes original movie
            
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
            obj.counter = 1;
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

        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages*obj.totalRuns;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages*obj.totalRuns;
        end
    end
end