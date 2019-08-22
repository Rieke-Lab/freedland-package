classdef (Abstract) RFDiskArrayProtocol < edu.washington.riekelab.protocols.RiekeLabProtocol
    
    properties (Access = protected)
        waitingForHardwareToStart
        changeRuns = 0
        runNo = 1
    end
    
    methods (Abstract)
        p = createPresentation(obj);
    end
    
    methods
        
        function redefineSettings(obj)
            
            % Predefined settings
            possibleoverrideRadii = repmat([0 0.75 2 3],18,1);
            possibleimageNo = [5 5 7 70 5 5 7 70 5 5 5 7 70 5 5 5 5 5]';
            possibleoverrideCoordinate = repmat({'RF'},18,1);
            possiblexSliceFrequency = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 2 0 1]';
            possibleySliceFrequency = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 2 0 1]';
            possiblerotateSlices = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 45 0 0 0]';
            possibledisksIgnoreCut = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0;...
                2 3; 0 3; 0 0; 0 0; 0 0; 2 3; 2 3; 2 3; 0 0; 0 2];
            possiblemeanDisks = [2 3 0; 1 0 0; 1 0 0; 1 0 0; 1 2 0; 1 2 3; 1 2 3; 1 2 3;...
                1 2 3; 1 2 3; 1 2 3; 1 2 3; 1 2 3; 1 0 0; 1 0 0; 1 0 0; 1 3 0; 1 2 3];
            possiblebackgroundDisks = [1 0 0; 2 3 0; 2 3 0; 2 3 0; 0 0 3; 0 0 0; 0 0 0; 0 0 0;...
                0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 2 3 0; 2 3 0; 2 3 0; 2 0 0; 0 0 0];

            % Assign variables
            obj.overrideRadii = possibleoverrideRadii(obj.runNo,:);
            obj.imageNo = possibleimageNo(obj.runNo,:);
            obj.overrideCoordinate = possibleoverrideCoordinate{obj.runNo,1};
            obj.xSliceFrequency = possiblexSliceFrequency(obj.runNo,:);
            obj.ySliceFrequency = possibleySliceFrequency(obj.runNo,:);  
            obj.rotateSlices = possiblerotateSlices(obj.runNo,:);
            obj.disksIgnoreCut = possibledisksIgnoreCut(obj.runNo,:);
            obj.meanDisks = possiblemeanDisks(obj.runNo,:);
            obj.backgroundDisks = possiblebackgroundDisks(obj.runNo,:);

            obj.runNo = obj.runNo + 1;
            obj.changeRuns = 1;
        end
        
        function prepareEpoch(obj, epoch)
            prepareEpoch@edu.washington.riekelab.protocols.RiekeLabProtocol(obj, epoch);
            
            obj.waitingForHardwareToStart = true;
            epoch.shouldWaitForTrigger = true;
            
            frameMonitor = obj.rig.getDevices('Frame Monitor');
            if ~isempty(frameMonitor)
                epoch.addResponse(frameMonitor{1});
            end
        end
        
        function controllerDidStartHardware(obj)
            controllerDidStartHardware@edu.washington.riekelab.protocols.RiekeLabProtocol(obj);
            
            if obj.playSequence == true
                if mod(obj.numEpochsCompleted,obj.numberOfAverages) == 0 && obj.numEpochsCompleted > 0
                    temp1 = obj.numEpochsPrepared;
                    temp2 = obj.numEpochsCompleted;
                    temp3 = obj.numIntervalsPrepared;
                    temp4 = obj.numIntervalsCompleted;
                    
                    prepareRun(obj)
                    
                    obj.numEpochsPrepared = temp1;
                    obj.numEpochsCompleted = temp2;
                    obj.numIntervalsPrepared = temp3;
                    obj.numIntervalsCompleted = temp4;
                end
            end

            if obj.waitingForHardwareToStart && obj.changeRuns == 1
                obj.changeRuns = 0;
                obj.waitingForHardwareToStart = false;
                obj.rig.getDevice('Stage').play(obj.createPresentation());
            elseif obj.waitingForHardwareToStart && obj.changeRuns == 0 % same run
                obj.waitingForHardwareToStart = false;
                obj.rig.getDevice('Stage').replay();
            end
        end
        
        function tf = shouldContinuePreloadingEpochs(obj) %#ok<MANU>
            tf = false;
        end
        
        function tf = shouldWaitToContinuePreparingEpochs(obj)
            tf = obj.numEpochsPrepared > obj.numEpochsCompleted || obj.numIntervalsPrepared > obj.numIntervalsCompleted;
        end
        
        function completeRun(obj)
            completeRun@edu.washington.riekelab.protocols.RiekeLabProtocol(obj);
        end
        
        function [tf, msg] = isValid(obj)
            [tf, msg] = isValid@edu.washington.riekelab.protocols.RiekeLabProtocol(obj);
            if tf
                tf = ~isempty(obj.rig.getDevices('Stage'));
                msg = 'No stage';
            end
        end
        
    end
    
end

