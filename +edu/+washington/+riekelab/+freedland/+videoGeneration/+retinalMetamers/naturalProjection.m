% Function within retinalMetamers
% By J. Freedland, 2020
%
%%% Project a naturalistic retinal movie into a low-dimensional
%%% representation. 
%
% Input:    obj: see DEMO
%           filename: name of exported MPG. Set to 0 to prevent exporting
%
% Output:   stimulus: structure with fields:
%               raw:        naturalistic movie
%               diskArray:  low-dimensional representation
%               values:     each region's corresponding light value.
%
% Last edited: 03/20 JF
%%%

function filename_adj = naturalProjection(obj, filename, varargin)

    ip = inputParser();
    ip.addParameter('exportRawMovie', false);
    ip.parse(varargin{:});

    disp('Building basic trajectory...')

    % Pull base trajectories and image information.
    [path,img,obj.saccades] = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.pathDOVES(obj.imageNo, obj.observerNo);
                
    % Normalize image to monitor
    img = (img./max(max(img)));
    img = img.*255;
    obj.imageMatrix = uint8(img);
    
    % Mean light intensity for retinal adaptation.
    obj.backgroundIntensity = mean(img(:));

    % Eye movement patterns from DOVES database.
    DOVES_trajectory    = 1/200:1/200:(length(path.x)/200); % 200 Hz
    monitorTrajectory   = 1/obj.monitorFrameRate:1/obj.monitorFrameRate:(length(path.x)/200);
    obj.xTraj = interp1(DOVES_trajectory,path.x,monitorTrajectory);
    obj.yTraj = interp1(DOVES_trajectory,path.y,monitorTrajectory);
    
    % Convert monitor to units in DOVES database.
    obj.videoSize = floor(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.monitorSize,obj.micronsPerPixel,'PIX2VH')) + 1;

    % Calculate individual neuron's receptive field (RF).
    [RFFilter,obj.rfSizing] = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.calculateFilter(obj);
    
    % Convolve filter with trajectory
    [weightedTraj, stimulus.raw] = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.weightedTrajectory(obj, img, RFFilter);
    
    disp('Calculating low-dimensional projection...')

    % Calculate disks
    [stimulus.projection, stimulus.values, stimulus.masks] = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.linearEquivalency(obj, weightedTraj, RFFilter, stimulus.raw);
    
    if obj.smoothing == true
        disp('Smoothing movies...')
        stimulus.projection = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.applySmoothing(stimulus.projection,stimulus.masks);
    end
    
    disp('Exporting movies...')
    if ip.Results.exportRawMovie == true
        filename_adj = strcat('img',string(obj.imageNo),'_raw');
        edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.exportMovie(obj,stimulus.raw,filename_adj)
    end
    
    if ~isequal(filename,0)
        filename_adj = strcat(filename,'_projection');
        edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.exportMovie(obj,stimulus.projection,filename_adj)
    end
    
    disp('Complete.')
end