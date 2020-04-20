% Uses dampening to find good disk sizes for retinal projections and
% metamers.
classdef spatialProjection < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    
    properties
        % Stimulus timing
        preTime = 250 % ms
        stimTime = 5500 % ms
        tailTime = 250 % ms
        
        % RF field information
        rfSigmaCenter = 70; % (um) enter from difference of gaussians fit for overlaying receptive field.
        rfSigmaSurround = 170; % (um) enter from difference of gaussians fit for overlaying receptive field.
        diskRadii = [75 200 400]; % radius for each disk, in microns. A cut will automatically be placed around the edge of the monitor.
        
        % Options
        naturalImages = [5,81];
        randomizeOrder = true
        numberOfAverages = uint16(5) % number of epochs to queue
        onlineAnalysis = 'extracellular'
        amp % Output amplifier
    end

    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})   

        % Rig options
        micronsPerPixel
        monitorSize
        monitorFrameRate
        diskRadii_arcmin
        videoSize
        
        % Assigned parameters
        imageNo  
        observerNo 
        numberOfDisks
        diskRegions
        diskRegionUnits
        slices
        cuts
        meanDisks 
        backgroundDisks
        switchDisks
        naturalDisks
        numberOfMetamerMovies
        smoothing
        switchPref
        saccades
        fixations
        imageMatrix
        backgroundIntensity
        xTraj
        yTraj
        rfSizing
        radii
        theta
        switchTraj
        frameNumber
        imageNumber
        
        % Load settings
        filenames
        directory
        sequence
        counter
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
            
            % Pull variables
            obj.micronsPerPixel = obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel');
            obj.monitorSize = obj.rig.getDevice('Stage').getCanvasSize();
            obj.monitorSize = fliplr(obj.monitorSize); % Adjust to [height, width]
            obj.monitorFrameRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
            obj.videoSize = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.monitorSize,obj.micronsPerPixel,'PIX2VH');
            
            % Make identify disk radii
            obj.diskRadii_arcmin = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.diskRadii,obj.micronsPerPixel,'UM2VH'));
            obj.diskRadii_arcmin = [0 obj.diskRadii_arcmin...
                round(max(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.monitorSize,obj.micronsPerPixel,'PIX2VH'))/2)];
            obj.diskRegionUnits = 'arcmin';
            
            % Create and export movies
            projectionNames = generateProjections(obj);

            metamerNames = generateMetamers(obj,2);
            compiledFilenames = [projectionNames; metamerNames];
            
            rawMovieNames = cell(length(obj.naturalImages),1);
            for a = 1:length(obj.naturalImages)
                rawMovieNames{a,1} = strcat('img',obj.naturalImages(a),'_raw');
            end
            
            % Disperse raw movie plays between movies
            frequency = 5; % how many movies to play between each raw movie?
           
            obj.filenames = [];
            c = 1;
            for a = 1:size(compiledFilenames,1)
                obj.filenames = [obj.filenames {compiledFilenames{a,1}}];
                
                if mod(a,frequency) == 0
                    obj.filenames = [obj.filenames {rawMovieNames{c,1}}];
                    c = mod(c,2) + 1;
                end
            end
            
            obj.counter = 1;
            obj.sequence = repmat(1:size(obj.filenames,1),1,obj.numberOfAverages);
            
            if obj.randomizeOrder == true
                obj.sequence = obj.sequence(randperm(length(obj.sequence)));
            end
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            epoch.addParameter('movieName',obj.filenames{obj.sequence(obj.counter),1});
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            p = stage.core.Presentation(obj.totalStimulusTime * 1e-3);

            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            obj.directory = 'Documents\freedland-package\+edu\+washington\+riekelab\+freedland\+movies';
            
            selection = obj.sequence(obj.counter);
            f = obj.filenames{selection,1}; % filename for relevant movie

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
        
        function filenames = generateProjections(obj)
            % This function uses user defined settings to generate, export,
            % and play a video within Symphony.
            experimentNumbers = [1:7,10,11,14];
            
            % Universal settings
            obj.observerNo  = 1;
            obj.smoothing   = false;
            obj.diskRegionUnits = 'arcmin';

            % Varying settings
            d   = [1,2,3,...
                   3,3,...
                   2,2,...
                   1,2,4];
            
            dRad = cell(1,length(d));
            for a = 1:length(d)
                dRad{1,a} = obj.diskRadii_arcmin(1:d(a)+1);
            end
                
            nDisks = [{1},{[1 2]},{[1 3]},...
                            {1},{1},...
                            {1},{1},...
                            {0},{0},{0}]; 
            bDisks = [{0},{0},{2},...
                               {2},{2},...
                               {0},{0},...
                               {0},{0},{0}]; 
            sDisks  = [{0},{0},{0},...
                            {3},{3},...
                            {2},{2},...
                            {0},{0},{0}]; 
            sPref   = [{''},{''},{''},...
                            {'darkFirst'},{'brightFirst'},...
                            {'darkFirst'},{'brightFirst'},...
                            {''},{''},{''}];
            mDisks  = [{0},{0},{0},...
                       {0},{0},...
                       {0},{0},...
                       {1},{[1 2]},{[1 2 3 4]}];
            slice     = [{0},{0},{0},...
                       {0},{0},...
                       {0},{0},...
                       {8},{4},{4}];
            dCut = [{0},{0},{0},...
                       {0},{0},...
                       {0},{0},...
                       {1},{1},{[1 2]}];
                        
            % Images to produce
            filenames = cell(length(experimentNumbers) * length(obj.naturalImages),1);
            counter1 = 1;
            
            for a = 1:length(obj.naturalImages)
                obj.imageNo = obj.naturalImages(a);
                                
                for b = 1:length(experimentNumbers)
                    obj.numberOfDisks       = d(b);
                    obj.diskRegions         = dRad{1,b};
                    obj.backgroundDisks     = bDisks{1,b}; 
                    obj.slices              = slice{1,b};
                    obj.cuts                = dCut{1,b};
                    obj.meanDisks           = mDisks{1,b};
                    obj.backgroundDisks     = bDisks{1,b};
                    obj.naturalDisks        = nDisks{1,b};
                    obj.switchDisks         = sDisks{1,b}; 
                    obj.switchPref          = sPref{1,b};
                    
                    if b == 1
                        edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.naturalProjection(obj,...
                            0,'exportRawMovie',true);
                    end
                    
                    filename = strcat('img',mat2str(obj.imageNo),'_',mat2str(obj.rfSigmaCenter),...
                        '_',mat2str(obj.rfSigmaSurround),'_',mat2str(experimentNumbers(b)));

                    if obj.smoothing == true
                        filename = strcat(filename,'_blur_');
                    end
                    
                    outputFilename = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.naturalProjection(obj,filename);
                    
                    filenames{counter1,1} = outputFilename;
                    counter1 = counter1 + 1;
                end
            end           
        end
        
        function filenames = generateMetamers(obj,numberMovies)
            % This function uses user defined settings to generate, export,
            % and play a video within Symphony.
            experimentNumbers = [8:9,12,13,15];
            smoothedExperiments = [12,15];

            obj.observerNo          = 1;  % Individual human observer number, as defined in the DOVES database (1-19)
            obj.backgroundDisks     = 0; 
            obj.switchDisks         = 0; 
            obj.naturalDisks        = 0;
            obj.numberOfMetamerMovies = numberMovies;
            obj.diskRegionUnits = 'arcmin';

            % Varying settings
            d   = [1,1,2,3,4];

            dRad = cell(1,length(d));
            for a = 1:length(d)
                dRad{1,a} = obj.diskRadii_arcmin(1:d(a)+1);
            end

            slice  = [0,4,4,4,4];
            cut    = [{0},{1},...
                       {[1 2]},{[1 2]},...
                       {[1 2 3 4]}];
            mDisks  = [{1},{1},...
                          {[1 2]},{[1 2 3]},...
                          {[1 2 3 4]}];

            % Images to produce
            filenames = cell(length(d) * length(obj.naturalImages) + length(smoothedExperiments),1);
            counter1 = 1;

            for a = 1:length(obj.naturalImages)
                obj.imageNo = obj.naturalImages(a);

                for b = 3:length(d)

                    obj.numberOfDisks	= d(b);
                    obj.diskRegions     = dRad{1,b}; % A vector that bounds each disk region.
                    obj.slices          = slice(b);
                    obj.cuts            = cut{1,b};
                    obj.meanDisks       = mDisks{1,b};
                    obj.smoothing       = false;

                    filename = strcat('img',mat2str(obj.imageNo),'_',mat2str(obj.rfSigmaCenter),...
                        '_',mat2str(obj.rfSigmaSurround),'_',mat2str(experimentNumbers(b)));

                    if sum(experimentNumbers(b) == smoothedExperiments) > 0
                        for c = 1:2 % Add additional smoothing

                            if c == 1 % Add blur
                                obj.smoothing = true;
                                filename = strcat(filename,'_blur_');
                            else      % Compare without blur
                                obj.smoothing = false;
                                filename = strcat('img',mat2str(obj.imageNo),'_',mat2str(obj.rfSigmaCenter),...
                                    '_',mat2str(obj.rfSigmaSurround),'_',mat2str(experimentNumbers(b)));
                            end

                            outputFilename = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.metamerGeneration(obj,filename,...
                                'exportProjection',true);

                            filenames{counter1,1} = outputFilename;
                            counter1 = counter1 + 1;
                        end
                    else
                        outputFilename = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.metamerGeneration(obj,filename,...
                                'exportProjection',true);

                            filenames{counter1,1} = outputFilename;
                            counter1 = counter1 + 1;
                    end    
                end
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