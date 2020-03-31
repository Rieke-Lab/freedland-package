% Replace a natural movie with a variety of integrated disks.
% By J. Freedland, 2019.
classdef RFMetamer < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Basics
        totalStimulusTime = 6000; % in ms
        referenceImage = 81; % image to build metamer for.
        
        % Arrangement of disks
        modelNumber = [15,13,12,9,8]; % Pre-determined model number to show replacement over.
        
        % RF field information
        rfSigmaCenter = 30; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.
        addBlur = true;
        randomize = true; % whether to randomize movies shown
        numberOfDistinctMovies = 3; % number of distinct movies to compare
        
        % Additional parameters
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

            obj.showFigure('symphonyui.builtin.figures.ResponseFigure', obj.rig.getDevice(obj.amp));
            obj.showFigure('edu.washington.riekelab.freedland.figures.MeanResponseFigure',...
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,'splitEpoch',1);
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
                        
            obj.directory = 'Documents\freedland-package\+edu\+washington\+riekelab\+freedland\+movies';
            D = dir(obj.directory);
            
            obj.movieFilenames = [];
            for b = 1:length(obj.modelNumber)
                % Identify relevant filenames
                specificFile = strcat('rep',mat2str(obj.referenceImage),'_',mat2str(obj.rfSigmaCenter),...,
                    '_',mat2str(obj.rfSigmaSurround),'_',mat2str(obj.modelNumber(b)));
                replacementMovies = [];
                
                for a = 1:size(D,1)
                    A = D(a).name;

                    % Only select relevant videos
                    if sum(strfind(A,'ref') & strfind(A,mat2str(obj.referenceImage))) > 0 % Normal trajectory
                        obj.movieFilenames = [obj.movieFilenames;{A}];
                    end

                    if sum(strfind(A,specificFile)) > 0 % Replacement trajectory
                        replacementMovies = [replacementMovies;{A}];
                    end
                end
            
                % Only consider blurred/non-blurred movies
                A = strfind(replacementMovies,'blur');
                if obj.addBlur == false
                    B = cellfun('isempty',A);
                else
                    B = not(cellfun('isempty',A));
                end
                replacementMovies = replacementMovies(B,:);

                % Check sufficient number of movies have been pre-generated
                if obj.numberOfDistinctMovies > size(replacementMovies,1)
                    msg = strcat('Not enough distinct movies pre-generated. Only',' ',mat2str(size(replacementMovies,1)),' movies exist.');
                    error(msg)
                else
                    replacementMovies = replacementMovies(1:obj.numberOfDistinctMovies,:);
                end
                obj.movieFilenames = [obj.movieFilenames;replacementMovies];
            end

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
            tf = obj.numEpochsPrepared < obj.numberOfAverages*obj.numberOfDistinctMovies;
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages*obj.numberOfDistinctMovies;
        end
    end
end