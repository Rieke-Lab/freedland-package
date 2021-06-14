% Flash rotated naturalistic images.
% By J. Freedland, 2019.
classdef flashImages < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime     = 250 % in ms
        stimTime    = 250 % in ms
        tailTime    = 250 % in ms
        
        % Region of image to be shown
        maskRadius  = 1000;  % in um. Places mask over surrounding portion of image.
        region      = 'centerWithAnnulus'; % where to display image
        
        % Natural image information
        cellPolarity = 'on'; % on- or off-pathway
        rectificationBias   = 'centerSurround'; % flash images with rectification only in the center, surround, or both. Select "all" to sample all three. Select "ignore" to simply sample a variety of luminances. Select "centerSurround" to focus on distibuted center/surround images.
        backgroundIntensity = 0.168; % common luminance to hold images at.
        
        % Reduce image appropriately
        includeReducedImages    = false;    % include 16-D reduce image flashes
        rfSigmaCenter           = 50;       % receptive-field center (only needed if reducing images)
        rfSigmaSurround         = 160;      % receptive-field surround (only needed if reducing images)

        % Additional parameters
        randomize   = true;  % randomize order of flashed images
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        regionType = symphonyui.core.PropertyType('char', 'row', {'center-only', 'surround-only', 'full-field','centerWithAnnulus','fullFieldWithAnnulus'}) 
        cellPolarityType = symphonyui.core.PropertyType('char', 'row', {'on', 'off'}) 
        rectificationBiasType = symphonyui.core.PropertyType('char', 'row', {'center','surround','full-field','all','ignore','centerSurround'}) 
        imageNo
        frame
        imageDatabase
        tracker
        counter
        order
        imageTracker
        frameTracker
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
            
            % Select dataset
            if strcmp(obj.cellPolarity,'on')
                if strcmp(obj.rectificationBias,'center')
                    % Surrounds entirely within [0, 50]% contrast, centers have >1
                    % regions outside [0, 50]%
                    obj.imageNo = [6,20,22,31,39,45,51,58,67,74,76,80,82,83,86,88,91];
                    obj.frame   = [289,574,289,635,64,281,127,1013,457,410,554,846,498,514,792,314,355];
                elseif strcmp(obj.rectificationBias,'surround')
                    % Centers entirely within [0, 50]% contrast, surrounds have >2
                    % regions outside [0, 50]%
                    obj.imageNo = [11,12,25,28,29,32,34,36,41,44,46,49,52,55,57,65,66,79,80,82,85,86,88,89,98];
                    obj.frame   = [357,393,299,423,589,202,56,70,615,399,262,176,164,954,143,340,159,676,244,877,690,621,745,62,265];
                elseif strcmp(obj.rectificationBias,'full-field')
                    % Centers have four regions within and outside [0, 50]% contrast
                    % Surrounds have four regions within and outside [0, 50]% contrast
                    obj.imageNo = [1,8,11,14,17,18,20,28,31,39,41,45,49,56,57,59,67,69,76,81,85,86,89,93];
                    obj.frame   = [495,971,488,758,797,422,926,270,744,593,118,805,359,489,309,277,636,129,341,731,59,248,775,570];
                elseif strcmp(obj.rectificationBias,'all')
                    % All datasets
                    obj.imageNo = [6,20,22,31,39,45,51,58,67,74,76,80,82,83,86,88,91,11,12,25,28,29,32,34,36,41,44,46,49,52,55,57,65,66,79,80,82,85,86,88,89,98,1,8,11,14,17,18,20,28,31,39,41,45,49,56,57,59,67,69,76,81,85,86,89,93];
                    obj.frame   = [289,574,289,635,64,281,127,1013,457,410,554,846,498,514,792,314,355,357,393,299,423,589,202,56,70,615,399,262,176,164,954,143,340,159,676,244,877,690,621,745,62,265,495,971,488,758,797,422,926,270,744,593,118,805,359,489,309,277,636,129,341,731,59,248,775,570];
                elseif strcmp(obj.rectificationBias,'ignore')
                    obj.imageNo = [1,84,17,77,43,89,20,60,85,86,100,2,72,42,25,28,28,27,28,3,37,93,53,76,88,82,67,80,20,11,98,83,31,37,11,9,32,66,28,48,15,67,52,24,14,9,76,67,37,6,67,99,6,57,99,25,97,14,49,21,48,14,70,15,52,60,5,2,97,101,31,40,56,90,99,6,18,62,65,101,77,99,89,99,56,31,9,88,31,18,26,9,26,36,99,97,54,5];
                    obj.frame = [57,1014,267,292,373,448,795,226,836,203,486,825,592,479,766,379,666,471,462,135,456,434,259,308,495,234,268,991,926,749,123,167,857,270,293,947,578,795,920,503,442,39,959,546,412,650,793,102,732,60,751,732,451,704,614,203,893,286,51,997,693,617,118,632,224,743,251,957,760,568,208,64,943,144,340,603,797,62,175,530,341,61,222,568,550,985,426,796,1022,201,202,135,234,878,450,857,254,190];
                elseif strcmp(obj.rectificationBias,'centerSurround')
                    % Images with evenly-distributed diversity of center-
                    % and surround-luminances
                    obj.imageNo = [1,19,92,98,98,83,64,30,7,2,5,4,42,8,73,86,81,93,93,46,29,73,72,77,93,96,94,47,100,85,25,27,78,77,8,91,69,90,22,91,28,12,67,53,97,49,4,83,86,19,57,14,18,45,9,74,68,66,15,98,74,2];
                    obj.frame = [225,626,871,569,504,923,719,996,182,863,363,901,888,542,404,152,912,955,833,856,110,524,421,1016,75,825,632,131,165,59,716,147,495,743,873,415,129,733,345,940,12,225,268,963,168,848,957,167,689,917,65,671,859,230,696,287,833,374,318,163,828,167];
                end
            elseif strcmp(obj.cellPolarity,'off')
                if strcmp(obj.rectificationBias,'center')
                    % Surrounds entirely within [-50, 0]% contrast, centers have >2
                    % regions outside [-50, 0]%
                    obj.imageNo = [7,8,13,16,18,19,23,33,35,38,44,47,75,77,78,80,83,93,95,98];
                    obj.frame   = [294,399,475,53,620,398,463,1016,599,164,636,409,536,341,451,538,923,732,693,569];
                elseif strcmp(obj.rectificationBias,'surround')
                    % Centers entirely within [-50, 0]% contrast, surrounds have >3
                    % regions outside [-50, 0]%
                    obj.imageNo = [3,5,6,10,12,20,21,23,25,29,34,40,51,58,64,74,79,83,84,85,86,90,92,96,98];
                    obj.frame   = [424,416,352,426,650,235,857,116,570,828,189,712,360,134,719,462,565,431,119,624,109,217,824,391,813];
                elseif strcmp(obj.rectificationBias,'full-field')
                    % Centers have four regions within and outside [-50, 0]% contrast
                    % Surrounds have four regions within and outside [-50, 0]% contrast
                    obj.imageNo = [3,8,12,31,33,45,46,67,69,73,85,87,89,93,94,100];
                    obj.frame   = [475,971,491,744,727,805,545,636,335,906,59,901,775,570,220,707];
                elseif strcmp(obj.rectificationBias,'all')
                    obj.imageNo = [7,8,13,16,18,19,23,33,35,38,44,47,75,77,78,80,83,93,95,98,3,5,6,10,12,20,21,23,25,29,34,40,51,58,64,74,79,83,84,85,86,90,92,96,98,3,8,12,31,33,45,46,67,69,73,85,87,89,93,94,100];
                    obj.frame   = [294,399,475,53,620,398,463,1016,599,164,636,409,536,341,451,538,923,732,693,569,424,416,352,426,650,235,857,116,570,828,189,712,360,134,719,462,565,431,119,624,109,217,824,391,813,475,971,491,744,727,805,545,636,335,906,59,901,775,570,220,707];
                elseif strcmp(obj.rectificationBias,'ignore')
                    obj.imageNo = [5,90,22,25,6,40,40,52,52,28,62,41,22,83,45,28,37,45,18,19,31,22,90,51,70,55,101,6,51,19,60,15,53,40,60,1,6,74,69,49,66,74,7,8,33,47,30,42,27,90,44,56,87,61,92,95,18,4,21,1,83,1,96,33,3,1,56,17,29,19,71,19,56,71,98,19,15,2,71,100,69,69,5,69,15,60,56,71,60,69];
                    obj.frame   = [190,144,788,247,561,320,272,959,558,920,845,488,484,514,445,12,402,230,157,340,338,49,879,127,582,189,72,664,51,851,604,967,317,830,397,225,905,929,929,848,159,571,754,128,669,230,520,636,906,385,749,53,277,764,261,629,620,538,908,57,805,779,282,161,256,603,370,474,305,148,136,54,310,811,698,246,127,482,307,294,53,477,681,1001,832,549,889,655,487,668];
                elseif strcmp(obj.rectificationBias,'centerSurround')
                    obj.imageNo = [1,19,92,98,98,83,64,30,7,2,5,4,42,8,73,86,81,93,93,46,29,73,72,77,93,96,94,47,100,85,25,27,78,77,8,91,69,90,22,91,28,12,67,53,97,49,4,83,86,19,57,14,18,45,9,74,68,66,15,98,74,2];
                    obj.frame = [225,626,871,569,504,923,719,996,182,863,363,901,888,542,404,152,912,955,833,856,110,524,421,1016,75,825,632,131,165,59,716,147,495,743,873,415,129,733,345,940,12,225,268,963,168,848,957,167,689,917,65,671,859,230,696,287,833,374,318,163,828,167];
                end
            end

            % Use DOVES database to identify images along trajectory
            imageFrame = findImage(obj);
            obj.imageDatabase = uint8(imageFrame);
            obj.tracker = repmat({'image'},size(imageFrame,4),1);
            obj.imageTracker = obj.imageNo;
            obj.frameTracker = obj.frame;

            if obj.includeReducedImages == true
                % Load generic settings
                settings = edu.washington.riekelab.freedland.videoGeneration.demoUtils.loadSettings(...
                    obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                    obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'),...
                    obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'),...
                    obj.rfSigmaCenter,obj.rfSigmaSurround);
                
                % Specific settings for reducing image
                settings.diskRegions = settings.diskRadii([1 3 5]); % Place one disk between regions [1] and [3]
                settings.meanDisks   = [1 2];                    % Make disk a linear-equivalent projection
                settings.slices      = 8;
                settings.sliceDisks  = [1 2];                    % Apply slicing to our single disk
                settings.backgroundIntensity = obj.backgroundIntensity*255;
                
                RFFilter = edu.washington.riekelab.freedland.videoGeneration.rfUtils.calculateFilter(settings);
                convolvedImage = imageFrame .* repmat(RFFilter,1,1,1,size(imageFrame,4));
                projection = edu.washington.riekelab.freedland.videoGeneration.utils.linearEquivalency(settings, convolvedImage, RFFilter, imageFrame);
                obj.imageDatabase = cat(4,obj.imageDatabase,uint8(projection));
                obj.tracker = cat(1,obj.tracker,repmat({'reduced'},size(projection,4),1));
                obj.imageTracker = cat(2,obj.imageNo,obj.imageNo);
                obj.frameTracker = cat(2,obj.frame,obj.frame);
            end

            % Setup display
            obj.counter = 0;
            if obj.randomize == true
                obj.order = randperm(size(obj.imageDatabase,4));
            else
                obj.order = 1:size(obj.imageDatabase,4);
            end
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('imageDisplayed',obj.imageTracker(obj.order(obj.counter+1)));
            epoch.addParameter('frameDisplayed',obj.frameTracker(obj.order(obj.counter+1)));
            epoch.addParameter('flashType',obj.tracker{obj.order(obj.counter+1),1});
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();     
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            
            % Rotate image
            specificImage = obj.imageDatabase(:,:,1,obj.order(obj.counter+1));
            p.setBackgroundColor(obj.backgroundIntensity)   % Set background intensity
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(specificImage);
            scene.size = [canvasSize(1),canvasSize(2)];
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);

            % Only allow image to be visible during specific time
            p.addStimulus(scene);
            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            % Create aperature
            aperatureDiameter = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.maskRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix') .* 2);
            aperture = stage.builtin.stimuli.Rectangle();
            aperture.position = canvasSize/2;
            aperture.color = obj.backgroundIntensity;
            aperture.size = [max(canvasSize) max(canvasSize)];
            
            % For annulus settings
            annulusSize = 50; % in microns
            annulusSize = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                annulusSize,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix'));
            
            if strcmp(obj.region,'center-only') %% Create aperture
                mask = stage.core.Mask.createCircularAperture(aperatureDiameter/max(canvasSize), 1024);
            elseif strcmp(obj.region,'centerWithAnnulus')
                mask = stage.core.Mask.createCircularAperture((aperatureDiameter-annulusSize)/max(canvasSize), 1024);
            elseif strcmp(obj.region,'surround-only')
                mask = stage.core.Mask.createAnnulus(0,aperatureDiameter/max(canvasSize), 1024);
            elseif strcmp(obj.region,'fullFieldWithAnnulus')
                mask = stage.core.Mask.createAnnulus((aperatureDiameter-annulusSize)/max(canvasSize),...
                    (aperatureDiameter+annulusSize)/max(canvasSize), 1024);
            elseif strcmp(obj.region,'full-field')
                mask = stage.core.Mask.createCircularAperture(1, 1024);
            end
            aperture.setMask(mask);
            p.addStimulus(aperture);
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end
        
        function [imageFrame,backgroundIntensities] = findImage(obj)
            % Relevant frame size
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize() / 2; % radius  
            pixelRange = round(edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(canvasSize,...
                obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'pix2arcmin'));

            imageFrame = zeros(pixelRange(2)*2+1,pixelRange(1)*2+1,1,length(obj.imageNo));
            backgroundIntensities = zeros(length(obj.imageNo),1);
            for a = 1:length(obj.imageNo)
                [path, image] = edu.washington.riekelab.freedland.videoGeneration.utils.pathDOVES(obj.imageNo(a), 1);
                image = image./max(image(:));
                backgroundIntensities(a) = mean(image(:));
                image = image.*255;      

                % Pull image
                imageFrame(:,:,1,a) = image(path.y(obj.frame(a))-pixelRange(2):path.y(obj.frame(a))+pixelRange(2),...
                    path.x(obj.frame(a))-pixelRange(1):path.x(obj.frame(a))+pixelRange(1));
            end
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.order);
        end
    end
end