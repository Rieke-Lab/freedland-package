% Flashes white noise to identify subunits.
% Methodology courtesy of "Nishal P. Shah, N. Brackbill, C. Rhoades, A. Tikidji-Hamburyan, G. Goetz, A. Litke, A. Sher, E. P. Simoncelli, E.J. Chichilnisky. Inference of Nonlinear Receptive Field Subunits with Spike-Triggered Clustering. ELife, 2020"
% By J. Freedland, 2020.
classdef flashWhiteNoise < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime  = 250 % in ms
        stimTime = 250 % in ms
        tailTime = 250 % in ms
        
        % White noise information
        centerRadius = 100;  % in um (region where white noise is displayed)
        pixelSize    = 20.8; % in um (Shah et. al. uses 20.8 & 41.6 um)
        uniqueFrames = 500;  % number of unique frames of white noise
        
        % Luminance information
        backgroundIntensity = 0.168; % light intensity of background; 0-1
        contrast            = 0.9;   % contrast for non-background regions; 0-1 
        
        % Additional parameters
        onlineAnalysis      = 'extracellular'
        numberOfAverages    = uint16(1) % number of epochs to queue
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        imageData
        counter
        order
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
            
            % Identify number of pixels to generate
            canvasSize = fliplr(obj.rig.getDevice('Stage').getCanvasSize());
            pix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.pixelSize,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            noiseSize = ceil(canvasSize/pix); % Number of unique noisy pixels
            
            % Randomly generate noise
            obj.imageData = rand([noiseSize,obj.uniqueFrames]); % For analysis
            
            % Rectify (evenly distributed between 0 to 1)
            A = (obj.imageData < 1/3);
            B = (obj.imageData > 1/3 & obj.imageData <= 2/3);
            C = (obj.imageData > 2/3);
            obj.imageData(A) = obj.backgroundIntensity .* (1-obj.contrast) .* 255; % Negative contrast
            obj.imageData(B) = obj.backgroundIntensity .* 255;                     % No change
            obj.imageData(C) = obj.backgroundIntensity .* (1+obj.contrast) .* 255; % Positive contrast
            
            % Only consider within associated region
            [xx,yy] = meshgrid(1:noiseSize(2),1:noiseSize(1));
            r = sqrt((xx - noiseSize(2)/2).^2 + (yy - noiseSize(1)/2).^2);
            centerRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.centerRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            cutoff = (r > ceil(centerRadiusPix/pix));
            
            % Restrict result
            obj.imageData(repmat(cutoff,1,1,size(obj.imageData,3))) = obj.backgroundIntensity .* 255;
            
            obj.counter = 0;
            obj.order = 1:size(obj.imageData,3);
            
            % Display est. stimulus time
            totalTime = (obj.uniqueFrames * obj.numberOfAverages) * ((obj.preTime + obj.stimTime + obj.tailTime + 1500)/1000);
            disp(strcat('approx stimulus time (+1.5 sec delay):',mat2str(round(totalTime/60)),' minutes'));
        end
            
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            epoch.addParameter('imageMatrix',obj.imageData(:,:,obj.order(obj.counter+1)));
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();     
            p = stage.core.Presentation((obj.preTime + obj.stimTime*60*1000 + obj.tailTime) * 1e-3);

            % Set image
            p.setBackgroundColor(obj.backgroundIntensity) % Set background intensity
            imageMatrix = obj.imageData(:,:,obj.order(obj.counter+1));
            
            % Prep to display image
            scene = stage.builtin.stimuli.Image(imageMatrix);
            scene.size = canvasSize;
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

            obj.counter = mod(obj.counter + 1,length(obj.order));
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.order);
        end
    end
end