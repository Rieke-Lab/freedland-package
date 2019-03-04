% adapted from Max Turner (2016).
% modified by JMF 2019.
classdef RFContrastReversingGrating < edu.washington.riekelab.protocols.RiekeLabStageProtocol

    properties
        preTime = 250 % ms
        stimTime = 2000 % ms
        tailTime = 250 % ms
        contrast = 0.9 % relative to mean (0-1)
        temporalFrequency = 4 % Hz
        radii = [0 100 200 300 400]; % input from RFDiskArray, in px
        barWidth = [5 10 20 40 80 160] % um
        rotation = 0; % deg
        backgroundIntensity = 0.5 % (0-1)
        randomizeOrder = false;
        onlineAnalysis = 'none'
        numberOfAverages = uint16(5) % number of epochs to queue
        amp
    end
    
    properties (Hidden)
        ampType
        onlineAnalysisType = symphonyui.core.PropertyType('char', 'row', {'none', 'extracellular', 'exc', 'inh'})
        barWidthSequence
        currentBarWidth
        radiiCounter
        outerRadii
        innerRadii
        barCounter
    end
       
    properties (Hidden, Transient)
        analysisFigure
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
                obj.rig.getDevice(obj.amp),'recordingType',obj.onlineAnalysis,...
                'groupBy',{'outerRadius'});
            obj.showFigure('edu.washington.riekelab.freedland.figures.FrameTimingFigure',...
                obj.rig.getDevice('Stage'), obj.rig.getDevice('Frame Monitor'));
            
            % Create bar width sequence.
            obj.barWidthSequence = obj.barWidth;
            obj.radiiCounter = 0;
            obj.barCounter = 0;
            
            if obj.randomizeOrder == true
                obj.barWidthSequence = obj.barWidthSequence(randperm(size(obj.barWidthSequence,2)));
            end
            
            obj.outerRadii = obj.radii(1,obj.radiiCounter + 2);
            obj.innerRadii = obj.radii(1,obj.radiiCounter + 1);
            obj.currentBarWidth = obj.barWidthSequence(obj.barCounter + 1);
        end

        function p = createPresentation(obj)
            canvasSize = obj.rig.getDevice('Stage').getCanvasSize();
            %convert from microns to pixels...
            currentBarWidthPix = obj.rig.getDevice('Stage').um2pix(obj.currentBarWidth);
            
            p = stage.core.Presentation((obj.preTime + obj.stimTime + obj.tailTime) * 1e-3); %create presentation of specified duration
            p.setBackgroundColor(obj.backgroundIntensity); % Set background intensity
            
            % Create grating stimulus.
            grate = stage.builtin.stimuli.Grating('square'); %square wave grating
            grate.orientation = obj.rotation;
            grate.size = [canvasSize(1), canvasSize(2)];
            grate.position = canvasSize/2;
            grate.spatialFreq = 1/(2*currentBarWidthPix); %convert from bar width to spatial freq
            grate.color = 2*obj.backgroundIntensity;
            
            %calc to apply phase shift s.t. a contrast-reversing boundary
            %is in the center regardless of spatial frequency. Arbitrarily
            %say boundary should be positve to right and negative to left
            %crosses x axis from neg to pos every period from 0
            zeroCrossings = 0:(grate.spatialFreq^-1):grate.size(1); 
            offsets = zeroCrossings-grate.size(1)/2; %difference between each zero crossing and center of texture, pixels
            [shiftPix, ~] = min(offsets(offsets>0)); %positive shift in pixels
            phaseShift_rad = (shiftPix/(grate.spatialFreq^-1))*(2*pi); %phaseshift in radians
            phaseShift = 360*(phaseShift_rad)/(2*pi); %phaseshift in degrees
            grate.phase = phaseShift; %keep contrast reversing boundary in center
            p.addStimulus(grate);
            %make it contrast-reversing
            if (obj.temporalFrequency > 0) 
                grateContrast = stage.builtin.controllers.PropertyController(grate, 'contrast',...
                    @(state)getGrateContrast(obj, state.time - obj.preTime/1e3));
                p.addController(grateContrast); %add the controller
            end
            
            function c = getGrateContrast(obj, time)
                c = obj.contrast.*sin(2 * pi * obj.temporalFrequency * time);
            end
              
            if (obj.innerRadii > 0) % Create mask
                mask = stage.builtin.stimuli.Ellipse();
                mask.position = canvasSize/2;
                mask.color = obj.backgroundIntensity;
                mask.radiusX = obj.innerRadii;
                mask.radiusY = obj.innerRadii;
                p.addStimulus(mask); %add mask
            end
            
            if  (obj.outerRadii > 0) % Create aperture
                aperture = stage.builtin.stimuli.Rectangle();
                aperture.position = canvasSize/2;
                aperture.color = obj.backgroundIntensity;
                aperture.size = [canvasSize(1)*2, canvasSize(1)*2];
                mask = stage.core.Mask.createCircularAperture(obj.outerRadii / max(obj.radii) / 2, 1024); %circular aperture
                aperture.setMask(mask);
                p.addStimulus(aperture); %add aperture
            end
            
            % hide during pre & post
            grateVisible = stage.builtin.controllers.PropertyController(grate, 'visible', ...
                @(state)state.time >= obj.preTime * 1e-3 && state.time < (obj.preTime + obj.stimTime) * 1e-3);
            p.addController(grateVisible);
            
            obj.radiiCounter = mod(obj.radiiCounter+1,(length(obj.radii)-1));
            obj.outerRadii = obj.radii(1,obj.radiiCounter+2);
            obj.innerRadii = obj.radii(1,obj.radiiCounter+1);
            
            if obj.radiiCounter == 0;
                obj.barCounter = mod(obj.barCounter + 1,length(obj.barCounter)+1);
            end
            
            obj.currentBarWidth = obj.barWidthSequence(obj.barCounter+1);
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabStageProtocol(obj, epoch);
            
            device = obj.rig.getDevice(obj.amp);
            duration = (obj.preTime + obj.stimTime + obj.tailTime) / 1e3;
            epoch.addDirectCurrentStimulus(device, device.background, duration, obj.sampleRate);
            epoch.addResponse(device);
           
            % bar greater than 1/2 aperture size -> just split field grating.
            % Allows grating texture to be the size of the aperture and the
            % resulting stimulus is the same...

            epoch.addParameter('currentBarWidth', obj.currentBarWidth);
            epoch.addParameter('innerRadius', obj.innerRadii);
            epoch.addParameter('outerRadius', obj.outerRadii);
        end
 
        function tf = shouldContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared < obj.numberOfAverages * length(obj.radii) * length(obj.barWidth);
        end
        
        function tf = shouldContinueRun(obj)
            tf = obj.numEpochsCompleted < obj.numberOfAverages * length(obj.radii) * length(obj.barWidth);
        end
        
        
    end
    
end