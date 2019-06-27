% Replace a natural image with a metamer.
% By J. Freedland, 2019.
classdef metamerFlash < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime = 250   % in ms
        stimTime = 250  % in ms
        tailTime = 250  % in ms
        
        % Natural image trajectory
        referenceImageNo = 5;    % natural image number (1 to 101)
        referenceFrameNumber = 0;% specific frame in a eye movement trajectory. set to zero to be random.
        numberMetamers = 2;      % number of metamers to view (up to 100)
        numberAntiMetamers = 2;  % number of anti-metamers to view (up to 100)
        numberRandom = 1;        % number of random images to view
        randomizeTrials = true;  % whether to randomize presentations
        metamerPercentile = 100; % sample from what percentile of matching images (100 = best, 1 = worst)
        antiMetamerPercentile = 80; % sample from what percentile of anti-matching images (100 = best anti-metamer, 1 = worst anti-metamer)
        
        
        % Additional parameters
        onlineAnalysis = 'extracellular'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        backgroundIntensity
        imageMatrices
        trialOrder
        counter
        typeTracker
        infoTracker
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
            
            % load metamer information
            filename = '+edu/+washington/+riekelab/+freedland/+images/metamerData.mat';
            load(filename)
            
            % adjust reference frame if needed.
            if obj.referenceFrameNumber == 0
                temp = randperm(1017);
                obj.referenceFrameNumber = temp(1);
            end
            
            % from analysis file, pull relevant frame and image information
            A = totalBest{obj.referenceFrameNumber, obj.referenceImageNo};
            B = totalWorst{obj.referenceFrameNumber, obj.referenceImageNo};
            
            % ensure trajectories are equivalent
            sensitiveA = A(:,3);
            sensitiveB = B(:,3);
            A(sensitiveA > 1017,:) = [];
            B(sensitiveB > 1017,:) = [];
            metamerQuartile = size(A,1) - round(size(A,1)*obj.metamerPercentile/100); % take from the 70th quartile, adds variability to sampling
            antiMetamerQuartile = round(size(B,1)*obj.antiMetamerPercentile/100); % take from the 70th quartile, adds variability to sampling

            % arrange
            metamerInfo = A(metamerQuartile+1:metamerQuartile+obj.numberMetamers,2:3); % sorted from best to least best
            antiMetamerInfo = B(antiMetamerQuartile-(obj.numberAntiMetamers-1):antiMetamerQuartile,2:3); % sorted from worst to least worst
            a = randperm(101); % images
            b = randperm(1000); % frames
            randomImgInfo = [a(1:obj.numberRandom)' b(1:obj.numberRandom)'];
            obj.backgroundIntensity = 0;
            
            % produce alternative images
            [obj.imageMatrices, obj.infoTracker, obj.typeTracker] = produceImage(obj,metamerInfo,antiMetamerInfo,randomImgInfo);

            if obj.randomizeTrials == true
                obj.trialOrder = randperm(size(obj.imageMatrices,3));
            else
                obj.trialOrder = 1:size(obj.imageMatrices,3);
            end
                
            obj.counter = 1;  
        end
        
        function prepareEpoch(obj, epoch)
            
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;

            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('backgroundIntensity', obj.backgroundIntensity);
            
            row = obj.trialOrder(obj.counter);
            epoch.addParameter('imageNo', obj.infoTracker(row,1));
            epoch.addParameter('frameNo', obj.infoTracker(row,2));
            epoch.addParameter('imageType', obj.typeTracker{row,1});
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset')); % in pixels
        end
        
        function p = createPresentation(obj)
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();             
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3);
            individualStim = obj.trialOrder(obj.counter);

            % Set background intensity
            p.setBackgroundColor(obj.backgroundIntensity);
            
            % Prep to display image
            imageMatrix = uint8(obj.imageMatrices(:,:,individualStim));
            scene = stage.builtin.stimuli.Image(imageMatrix);
            
            % Use linear interpolation when scaling the image
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            scene.size = [size(imageMatrix,2) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),...
                size(imageMatrix,1) * 3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')];
            p0 = canvasSize/2;
            scene.position = p0;
            
            p.addStimulus(scene);

            sceneVisible = stage.builtin.controllers.PropertyController(scene, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(sceneVisible);
            
            obj.counter = mod(obj.counter,length(obj.trialOrder)) + 1;
        end
        
        function [outputImages, imageInfo, types] = produceImage(obj,metamers,antiMetamers,randomImgInfo)
            imageNumbs = [obj.referenceImageNo; metamers(:,1); antiMetamers(:,1); randomImgInfo(:,1)];
            frameNumbs = [obj.referenceFrameNumber; metamers(:,2); antiMetamers(:,2); randomImgInfo(:,2)];
            imageInfo = [imageNumbs frameNumbs];
            types = cell(size(imageNumbs));
            types{1} = 'original';
            
            for a = 2:size(metamers,1)+1
                types{a} = 'metamer';
            end
            
            baseline = size(metamers,1)+2;
            
            for a = baseline:baseline+size(antiMetamers,1)
                types{a} = 'anti-metamer';
            end
            
            baseline = size(metamers,1)+size(antiMetamers,1)+2;
            
            for a = baseline:baseline+size(randomImgInfo,1)
                types{a} = 'random';
            end
            
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            imgSize = ceil(canvasSize / (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel')));
            xRange = floor(imgSize(1) / 2);
            yRange = floor(imgSize(2) / 2);
            centering = obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'); % in mu
            centeringPix = centering ./ (3.3/obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            
            outputImages = zeros(imgSize(2),imgSize(1),size(imageNumbs,1));
 
            for a = 1:size(imageNumbs,1)
                imageVal = imageNumbs(a,1);
                frameVal = frameNumbs(a,1);
                [~, baseMovement, fixMovement, pictureInformation] = edu.washington.riekelab.freedland.scripts.pathDOVES(imageVal, 1,...
                        'amplification', 1,'mirroring', true); 
                img = pictureInformation.image;
                img = img./max(img(:)) .* 255;
                
                if a == 1
                    obj.backgroundIntensity = mean(img(:)) / 255;
                end
                
                xTraj = baseMovement.x(frameVal) + fixMovement.x(frameVal);
                yTraj = baseMovement.y(frameVal) + fixMovement.y(frameVal);

                % while the rig automatically centers the stimulus, our
                % calculation doesn't.
                centeredXTraj = round(xTraj - centeringPix(1));
                centeredYTraj = round(yTraj + centeringPix(2));

                % create image
                outputImages(:,:,a) = img(centeredYTraj-yRange:centeredYTraj+yRange-1,...
                    centeredXTraj-xRange:centeredXTraj+xRange-1); 
            end
        end 
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.trialOrder);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.trialOrder);
        end
    end
end