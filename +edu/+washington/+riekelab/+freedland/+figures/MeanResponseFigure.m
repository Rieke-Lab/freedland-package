classdef MeanResponseFigure < symphonyui.core.FigureHandler
    
    properties (SetAccess = private)
        device
        groupBy
        sweepColor
        recordingType
        storedSweepColor
        splitEpoch
    end
    
    properties (Access = private)
        axesHandle
        sweeps
        sweepIndex
        storedSweep
    end
    
    methods
        
        function obj = MeanResponseFigure(device, varargin)
            co = get(groot, 'defaultAxesColorOrder');

            ip = inputParser();
            ip.addParameter('groupBy', [], @(x)iscellstr(x));
            ip.addParameter('sweepColor', co(1,:), @(x)ischar(x) || ismatrix(x));
            ip.addParameter('storedSweepColor', 'r', @(x)ischar(x) || isvector(x));
            ip.addParameter('recordingType', [], @(x)ischar(x));
            ip.addParameter('splitEpoch', 1);
            ip.parse(varargin{:});
            
            obj.device = device;
            obj.groupBy = ip.Results.groupBy;
            obj.sweepColor = ip.Results.sweepColor;
            obj.storedSweepColor = ip.Results.storedSweepColor;
            obj.recordingType = ip.Results.recordingType;
            obj.splitEpoch = ip.Results.splitEpoch;
            
            obj.createUi();
        end
        
        function createUi(obj)
            import appbox.*;
            iconDir = [fileparts(fileparts(mfilename('fullpath'))), '\+utils\+icons\'];
            toolbar = findall(obj.figureHandle, 'Type', 'uitoolbar');
            
            pVEButton = uipushtool( ...
                'Parent', toolbar, ...
                'TooltipString', 'Percent variance explained', ...
                'Separator', 'on', ...
                'ClickedCallback', @obj.onSelectedpVE);
            setIconImage(pVEButton, [iconDir, 'DoG.png']);
            
            storeSweepButton = uipushtool( ...
                'Parent', toolbar, ...
                'TooltipString', 'Store Sweep', ...
                'Separator', 'on', ...
                'ClickedCallback', @obj.onSelectedStoreSweep);
            setIconImage(storeSweepButton, symphonyui.app.App.getResource('icons/sweep_store.png'));
            
            clearStoredButton = uipushtool( ...
                'Parent', toolbar, ...
                'TooltipString', 'Clear saved sweep', ...
                'Separator', 'off', ...
                'ClickedCallback', @obj.onSelectedClearStored);
            setIconImage(clearStoredButton, [iconDir, 'Xout.png']);
            
            obj.axesHandle = axes( ...
                'Parent', obj.figureHandle, ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            xlabel(obj.axesHandle, 'sec');
            obj.sweeps = {};
            obj.setTitle([obj.device.name ' Mean Response']);
        end
        
        function setTitle(obj, t)
            set(obj.figureHandle, 'Name', t);
            title(obj.axesHandle, t);
        end
        
        function clear(obj)
            cla(obj.axesHandle);
            obj.sweeps = {};
        end
        
        function handleEpoch(obj, epoch)
            if ~epoch.hasResponse(obj.device)
                error(['Epoch does not contain a response for ' obj.device.name]);
            end
            
            response = epoch.getResponse(obj.device);
            [quantities, units] = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            if numel(quantities) > 0
                x = (1:numel(quantities)) / sampleRate;
                y = quantities;
                
                if strcmp(obj.recordingType,'extracellular')
                    filterSigma = (15/1000)*sampleRate; %15 msec -> dataPts
                    newFilt = normpdf(1:10*filterSigma,10*filterSigma/2,filterSigma);
                    res = edu.washington.riekelab.freedland.utils.spikeDetectorOnline(y,[],sampleRate);
                    y = zeros(size(y));
                    y(res.sp) = 1; %spike binary
                    y = sampleRate*conv(y,newFilt,'same'); %inst firing rate, Hz
                end    
                
                g = floor(length(y)/obj.splitEpoch);
                y1 = y(1:g);
                
                if obj.splitEpoch > 1
                    y2 = y(g+1:2*g);
                end
                
                if obj.splitEpoch > 2
                    y3 = y(2*g+1:3*g);
                end
                
                x = x(1:g);
            else
                x = [];
                y = [];
            end
            
            p = epoch.parameters;
            if isempty(obj.groupBy) && isnumeric(obj.groupBy)
                parameters = p;
            else
                parameters = containers.Map();
                for i = 1:length(obj.groupBy)
                    key = obj.groupBy{i};
                    parameters(key) = p(key);
                end
            end
            
            if isempty(parameters)
                t = 'All epochs grouped together';
            else
                t = ['Grouped by ' strjoin(parameters.keys, ', ')];
            end
            obj.setTitle([obj.device.name ' Mean Response (' t ')']);
            
            obj.sweepIndex = [];
            for i = 1:numel(obj.sweeps)
                if isequal(obj.sweeps{i}.parameters, parameters)
                    obj.sweepIndex = i;
                    break;
                end
            end
            
            if isempty(obj.sweepIndex)
                if size(obj.sweepColor,1) == 1
                    cInd = 1;
                elseif size(obj.sweepColor,1) >= length(obj.sweeps)+1
                    cInd = length(obj.sweeps)+1;
                else
                    cInd = 1;
                    warning('Not enough colors supplied for sweeps')
                end
                
                if obj.splitEpoch == 2
                    sweep.line = line(x, y1, 'Parent', obj.axesHandle,...
                        'Color', 'blue');
                    sweep.line2 = line(x, y2, 'Parent', obj.axesHandle,...
                        'Color', 'red');
                elseif obj.splitEpoch == 3
                    sweep.line = line(x, y1, 'Parent', obj.axesHandle,...
                        'Color', 'blue');
                    sweep.line2 = line(x, y2, 'Parent', obj.axesHandle,...
                        'Color', 'red');
                    sweep.line3 = line(x, y3, 'Parent', obj.axesHandle,...
                        'Color', 'green');
                else
                    sweep.line = line(x, y, 'Parent', obj.axesHandle,...
                        'Color', obj.sweepColor(cInd,:));
                end
                sweep.parameters = parameters;
                sweep.count = 1;
                obj.sweeps{end + 1} = sweep;
            else
                sweep = obj.sweeps{obj.sweepIndex};
                cy = get(sweep.line, 'YData');
                if obj.splitEpoch == 2
                    cy2 = get(sweep.line2, 'YData');
                    set(sweep.line, 'YData', (cy * sweep.count + y1) / (sweep.count + 1));
                    set(sweep.line2, 'YData', (cy2 * sweep.count + y2) / (sweep.count + 1));
                elseif  obj.splitEpoch == 3
                    cy2 = get(sweep.line2, 'YData');
                    cy3 = get(sweep.line3, 'YData');
                    set(sweep.line, 'YData', (cy * sweep.count + y1) / (sweep.count + 1));
                    set(sweep.line2, 'YData', (cy2 * sweep.count + y2) / (sweep.count + 1));
                    set(sweep.line3, 'YData', (cy3 * sweep.count + y3) / (sweep.count + 1));
                else
                    set(sweep.line, 'YData', (cy * sweep.count + y) / (sweep.count + 1));
                end
                sweep.count = sweep.count + 1;
                obj.sweeps{obj.sweepIndex} = sweep;
            end
            
            %check for stored data to plot...
            storedData = obj.storedAverages();
            if ~isempty(storedData)
                if ~isempty(obj.storedSweep) %Handle still there
                    if obj.storedSweep.line.isvalid %Line still there
                        
                    else
                        obj.storedSweep.line = line(storedData(1,:), storedData(2,:),...
                        'Parent', obj.axesHandle, 'Color', obj.storedSweepColor);
                    end                 
                else %no handle
                    obj.storedSweep.line = line(storedData(1,:), storedData(2,:),...
                        'Parent', obj.axesHandle, 'Color', obj.storedSweepColor);
                end
            end

            ylabel(obj.axesHandle, units, 'Interpreter', 'none');
        end
        
    end
    
    methods (Access = private)
        
        function onSelectedStoreSweep(obj, ~, ~)
            if isempty(obj.sweepIndex)
                sweepPull = 1;
            else
                sweepPull = obj.sweepIndex;
            end
            if ~isempty(obj.storedSweep) %Handle still there
                if obj.storedSweep.line.isvalid %Line still there
                    %delete the old storedSweep
                    obj.onSelectedClearStored(obj)
                end
            end
            
            %save out stored data
            obj.storedSweep.line = obj.sweeps{sweepPull}.line;
            obj.storedAverages([obj.storedSweep.line.XData; obj.storedSweep.line.YData]);
            %set the saved trace to storedSweepColor to indicate that it has been saved
            obj.storedSweep.line = line(obj.storedSweep.line.XData, obj.storedSweep.line.YData,...
                        'Parent', obj.axesHandle, 'Color', obj.storedSweepColor);
        end

        function onSelectedClearStored(obj, ~, ~)
            obj.storedAverages('Clear');
            obj.storedSweep.line.delete
        end

    end
    
    methods (Static)
        
        
        function averages = storedAverages(averages)
            % This method stores means across figure handlers.
            persistent stored;
            if (nargin == 0) %retrieve stored data
               averages = stored;
            else %set or clear stored data
                if strcmp(averages,'Clear')
                    stored = [];
                else
                    stored = averages;
                    averages = stored;
                end
            end
        end
    end
    
    methods (Access = private)
        
        function onSelectedpVE(obj, ~, ~)
            if obj.splitEpoch == 2 % must have reference

                sigmaC = var(obj.sweeps{1,1}.line.YData - obj.sweeps{1,1}.line2.YData) / var(obj.sweeps{1,1}.line.YData);

                if sigmaC > 1
                    sigmaC = 1; % prevent from producing negative value.
                end
                
                str = {['percent variance explained = ', num2str((1-sigmaC)*100),'%']};
                title(obj.axesHandle,str);
            end
        end

    end
        
end

