% Offsets representations to measure error.
% By J. Freedland, 2019.
classdef offsetVideos < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % RF field information
        rfSigmaCenter = 30; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 100; % (um) enter from difference of gaussians fit for overlaying receptive field.

        imageNumber = 5;
        controlVideo = '1A'; % partial string contained in movie file
        testVideo = '1S-9'; %  partial string contained in movie file
        offsets = [0 25 50 75 100 150 200 400]; % offsets (um)
        randomize = true; % whether to randomize movies shown

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
        offsetsPix
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

            % Find movie names
            controlMov = [];
            testMov = [];
            imgStr = strcat('img',num2str(obj.imageNumber));
            for a = 1:length(D)
                if sum(strfind(D(a).name,obj.controlVideo)) > 0 || sum(strfind(D(a).name,imgStr)) > 0
                    controlMov = D(a).name;
                elseif sum(strfind(D(a).name,obj.testVideo)) > 0 || sum(strfind(D(a).name,imgStr)) > 0
                    testMov = D(a).name;
                end
            end
            
            if isempty(controlMov) || isempty(testMov)
                error('Cannot find correct movie.')
            end
            obj.movieFilenames = repmat([controlMov;testMov],length(obj.offsets),1);

            % Find background intensity
            [~, ~, ~, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(81, 1,...
                    'amplification', 1,'mirroring', false);
            img = pictureInformation.image;
            img = (img./max(max(img)));
            obj.backgroundIntensity = mean(img(:));
            
            obj.offsetsPix = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(...
                obj.offsets,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'UM2PIX');
            obj.offsetsPix = repelem(obj.offsetsPix,1,2);
            obj.sequence = repmat(1:length(obj.offsetsPix),1,obj.numberOfAverages);
            
            if obj.randomize == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
            obj.counter = 1;
        
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = 6;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('movieName',obj.movieFilenames{obj.sequence(obj.counter),1});
            
            offsetUm = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.changeUnits(...
                obj.offsetsPix(obj.sequence(obj.counter)),obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'UM2PIX');
            epoch.addParameter('offset',offsetUm);

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
            p0 = p0(1) + obj.offsetsPix(obj.counter);
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