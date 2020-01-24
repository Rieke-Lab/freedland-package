classdef (Abstract) RFDiskArrayProtocol < edu.washington.riekelab.protocols.RiekeLabProtocol
    
    properties (Access = protected)
        waitingForHardwareToStart
        changeRuns = 0
        runNo = 1
        firstRun = 1
        order
        numberOfExperiments = 16; % How many stimuli to create and test (in total).
    end
    
    methods (Abstract)
        p = createPresentation(obj);
    end
    
    methods
        
        function redefineSettings(obj,t)
            
            % t: whether to advance forward to the next setting (t/f).
            if obj.firstRun == 1 
                if obj.randomize == true
                    obj.order = randperm(obj.numberOfExperiments);
                else
                    obj.order = 1:obj.numberOfExperiments;
                end
                obj.firstRun = 0;
            end
            
            % Each row represents a different stimulus.
            possibleImageNo = repelem(81,obj.numberOfExperiments,1);
            
            possibleDisks = repelem(3,1,obj.numberOfExperiments);
            possibleDisks(9) = 4;
            
            possibleOverrideRadii = cell(obj.numberOfExperiments,1);
            possibleOverrideRadii(:) = {[0 0.75 2 3]};
            possibleOverrideRadii(9) = {[0 0.75 2 2.2 3]};

            possibleXSliceFrequency = [0 0 0,...
                                        0 0,...
                                        0 0 0 0,...
                                        0 0 1 2,...
                                        0 1 1];
            possibleYSliceFrequency = [0 0 0,...
                                        0 0,...
                                        0 0 0 0,...
                                        0 1 1 2,...
                                        0 1 1];
            possibleDisksIgnoreCut = [0 0; 0 0; 0 0; 
                                      0 0; 0 0; ...
                                      0 0; 0 0; 0 0; 0 0;...
                                      2 3; 2 3; 2 3; 2 3;...
                                      0 0; 2 3; 0 0];
            possibleMeanDisks = [0 0 0; 0 0 0; 0 0 0; 
                                0 0 0; 0 0 0;...
                                 1 0 0; 1 2 0; 1 0 3; 1 3 4;...
                                 1 0 0; 1 0 0; 1 0 0; 1 0 0;...
                                 1 2 3; 1 2 3; 1 2 3];
            possibleBackgroundDisks = [0 2 3; 0 0 3; 0 2 0; 
                                        0 2 0; 0 2 0;...
                                       0 2 3; 0 0 3; 0 2 0; 0 2 0;...
                                       0 2 3; 0 2 3; 0 2 3; 0 2 3;...
                                       0 0 0; 0 0 0; 0 0 0];
            possibleNaturalDisks = [1 0 0; 1 2 0; 1 0 3; 
                                    1 0 0; 1 0 0;...
                                    0 0 0; 0 0 0; 0 0 0; 0 0 0;...
                                    0 0 0; 0 0 0; 0 0 0; 0 0 0;...
                                    0 0 0; 0 0 0; 0 0 0];
            possibleSwitchDisks = [0 0 0; 0 0 0; 0 0 0; ...
                                    0 0 3; 0 0 3;...
                                    0 0 0; 0 0 0; 0 0 0; 0 0 0;...
                                    0 0 0; 0 0 0; 0 0 0; 0 0 0;...
                                    0 0 0; 0 0 0; 0 0 0];
            % More information can be added from RFDiskArray as needed.

            % Identify and assign parameters based on experiment #.
            runVal = obj.order(obj.runNo);
                
            obj.imageNo             = possibleImageNo(runVal);
            obj.disks               = possibleDisks(runVal);
            obj.overrideRadii       = possibleOverrideRadii{runVal,1};
            obj.xSliceFrequency     = possibleXSliceFrequency(runVal);
            obj.ySliceFrequency     = possibleYSliceFrequency(runVal);  
            obj.disksIgnoreCut      = possibleDisksIgnoreCut(runVal,:);
            obj.meanDisks           = possibleMeanDisks(runVal,:);
            obj.backgroundDisks     = possibleBackgroundDisks(runVal,:);
            obj.naturalDisks        = possibleNaturalDisks(runVal,:);
            obj.switchDisks         = possibleSwitchDisks(runVal,:);

            if t == true
                disp(obj.runNo);
                obj.runNo = obj.runNo + 1;
                obj.changeRuns = 1;
            end
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

