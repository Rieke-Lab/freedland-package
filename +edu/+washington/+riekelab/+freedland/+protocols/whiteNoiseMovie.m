% Flashes white noise to identify subunits.
% Methodology courtesy of "Nishal P. Shah, N. Brackbill, C. Rhoades, A. Tikidji-Hamburyan, G. Goetz, A. Litke, A. Sher, E. P. Simoncelli, E.J. Chichilnisky. Inference of Nonlinear Receptive Field Subunits with Spike-Triggered Clustering. ELife, 2020"
% By J. Freedland, 2020.
classdef whiteNoiseMovie < edu.washington.riekelab.protocols.RiekeLabStageProtocol
    properties
        % Stimulus timing
        preTime     = 250  % in ms
        tailTime    = 250  % in ms
        stimTime    = 10   % individual movie length, in seconds
        totalTime   = 15  % total recording time across all movies, in minutes
        
        % White noise information
        centerRadius = 100;  % in um (region where white noise is displayed)
        pixelSize    = 20.8; % in um (Shah et. al. uses 20.8 & 41.6 um)
        
        % Luminance information
        backgroundIntensity = 0.168; % light intensity of background; 0-1
        contrast            = 0.9;   % contrast for non-background regions; 0-1 
        
        % Additional parameters
        onlineAnalysis      = 'extracellular'
        amp % Output amplifier
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'}) 
        movieName
        ranges
        directory
        counter
        order
        movieMatrix
        filename
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
            
            % Only consider pixels within associated region
            [xx,yy] = meshgrid(1:noiseSize(2),1:noiseSize(1));
            r = sqrt((xx - noiseSize(2)/2).^2 + (yy - noiseSize(1)/2).^2);
            centerRadiusPix = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(...
                obj.centerRadius,obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'),'um2pix');
            cutoff = (r > ceil(centerRadiusPix/pix));
            
            % Randomly generate noise
            refreshRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
                
            % Generate noise as N x M x 1 x F matrix
            noise = rand([noiseSize,1,obj.totalTime*60*refreshRate]);
                
            % Evenly rectify into pixels
            A = (noise < 1/3);
            B = (noise > 1/3 & noise <= 2/3);
            C = (noise > 2/3);
            noise(A) = obj.backgroundIntensity .* (1-obj.contrast) .* 255; % Negative contrast
            noise(B) = obj.backgroundIntensity .* 255;                     % No change
            noise(C) = obj.backgroundIntensity .* (1+obj.contrast) .* 255; % Positive contrast
            noise(repmat(cutoff,1,1,1,size(noise,4))) = obj.backgroundIntensity .* 255; % Surround is static

            % Export data into .mat file
            obj.filename = strcat('Documents/freedland-package/+edu/+washington/+riekelab/+freedland/+movies/whiteNoise');
            obj.movieMatrix = noise;
            obj.movieName = strcat(obj.filename,'.mp4');
            save(obj.filename,'noise')
            
            % Define ranges for individual movies
            obj.ranges = 1:obj.stimTime*refreshRate:obj.totalTime*60*refreshRate;
            obj.counter = 0;
            obj.order = 1:length(obj.ranges);
        end
            
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
            
            % Frames displayed
            epoch.addParameter('range',[obj.ranges(obj.order(obj.counter+1)), obj.ranges(obj.order(obj.counter+2))]);
            
            % Add metadata from Stage, makes analysis easier.
            epoch.addParameter('canvasSize',obj.rig.getDevice('Stage').getConfigurationSetting('canvasSize'));
            epoch.addParameter('micronsPerPixel',obj.rig.getDevice('Stage').getConfigurationSetting('micronsPerPixel'));
            epoch.addParameter('monitorRefreshRate',obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            epoch.addParameter('centerOffset',obj.rig.getDevice('Stage').getConfigurationSetting('centerOffset'));
        end
        
        function p = createPresentation(obj)
            
            % Stage presets
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();     
            p = stage.core.Presentation((obj.preTime + obj.tailTime) * 1e-3 + obj.stimTime);

            % Set image
            p.setBackgroundColor(obj.backgroundIntensity) % Set background intensity
            
            % Define frames to show
            individualMovie = obj.movieMatrix(:,:,1,...
                obj.ranges(obj.order(obj.counter+1)):obj.ranges(obj.order(obj.counter+2)));
            
            % Scale to monitor
            individualMovie = imresize(individualMovie,[canvasSize(2) canvasSize(1)],'nearest');
            
            % Add null frames
            preFrames = round(obj.preTime / 1000 * obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            postFrames = round(obj.tailTime / 1000 * obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate'));
            nullFrame = ones(size(individualMovie,1),size(individualMovie,2)) .* obj.backgroundIntensity .* 255;
            individualMovie = uint8(cat(4,repmat(nullFrame,1,1,1,preFrames),...
                individualMovie,repmat(nullFrame,1,1,1,postFrames)));
            
            % Export as .mp4 for stage.
            v = VideoWriter(obj.filename,'MPEG-4');
            v.FrameRate = obj.rig.getDevice('Stage').getConfigurationSetting('monitorRefreshRate');
            open(v)
            for b = 1:size(individualMovie,4)
                writeVideo(v,individualMovie(:,:,:,b))
            end
            
            % Prep to display movie
            scene = stage.builtin.stimuli.Movie('Documents/freedland-package/+edu/+washington/+riekelab/+freedland/+movies/whiteNoise.mp4');
            scene.size = canvasSize;
            p0 = canvasSize/2;
            scene.position = p0;
            
            % Use linear interpolation when scaling the image (for small
            % changes)
            scene.setMinFunction(GL.LINEAR);
            scene.setMagFunction(GL.LINEAR);
            
            p.addStimulus(scene);
            obj.counter = mod(obj.counter + 1,length(obj.order));
        end
        
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < length(obj.order);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < length(obj.order);
        end
    end
end