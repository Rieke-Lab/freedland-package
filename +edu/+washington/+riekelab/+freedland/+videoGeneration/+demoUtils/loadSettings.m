% Load generic settings.

function settings = loadSettings(umPerPixel,canvasSize,frameRate,centerSigma,surroundSigma)

    % Cell information
    settings.rfSigmaCenter       = centerSigma;   % in microns
    settings.rfSigmaSurround     = surroundSigma; % in microns
    
    % Pixel information
    settings.micronsPerPixel     = umPerPixel; % How many microns each pixel spans
    settings.monitorSize         = fliplr(canvasSize); % Size of monitor [height, width]
    settings.monitorFrameRate    = frameRate;  % Hz
    settings.videoSize = edu.washington.riekelab.freedland.videoGeneration.utils.changeUnits(settings.monitorSize,settings.micronsPerPixel,'pix2arcmin');
    if mod(settings.videoSize,2) == 0
        settings.videoSize = settings.videoSize + 1; % Ensure is odd (s.t. central pixel exists)
    end

    % Calculate the location of disks
    [~,info] = edu.washington.riekelab.freedland.videoGeneration.rfUtils.calculateFilter(settings);
    settings.diskRadii = [0,...                                                % At center of monitor
        info.percentExcitation(2,(info.percentExcitation(1,:) == 50)),... % At 50% excitation
        info.zeroPt,...                                                   % Where excitatory inputs -> inhibitory
        info.percentInhibition(2,(info.percentInhibition(1,:) == 50)),... % At 50% inhibition
        max(settings.videoSize/2)];                                            % At edge of monitor
    settings.diskRegionUnits = 'arcmin';

    % Stimulus timing
    settings.preTime     = 250;  % (in ms). Presents blank background.
    settings.stimTime    = 5500; % (in ms). Presents main stimulus.
    settings.tailTime    = 250;  % (in ms). Presents blank background.
    
    %%% All settings
    % DOVES information
    settings.imageNo     = 5; % Individual image to show. (#1 - 101)
    settings.observerNo  = 1; % Individual observer for eye tracking. (#1 - 19)
    
    % List of all key variables
    settings.experimentName      = 'test';  % ID for corresponding movie (string)
    settings.diskRegions         = []; % Specific radii for placing disks
    
    % Different types of disks (see EXAMPLES below)
    % Here, we identify each disk from center outwards as #s 1,2,3,...
    settings.meanDisks           = []; % Replace disk #s with linear equivalent disk.
    settings.backgroundDisks     = []; % Replace disk #s with static disk (average intensity of image) 
    settings.naturalDisks        = []; % Replace disk #s with original image
    settings.switchDisks         = []; % Replace disk #s with flashing region (after each saccade)
        settings.switchContrast  = []; % Intensity of flashing region (0 - 1)
    settings.metamerDisks        = []; % Replace disk #s with another naturalistic image.

    %%% EXAMPLES
    %   - Make the center-most disk a linear equivalent disk: 
    %       obj.meanDisks = 1;
    %   - Make the two center-most disks a natural image: 
    %       obj.naturalDisks = [1 2];
    %   - Make the central disk a metamer and second disk switch brightness periodically.
    %       obj.metamerDisks = 1;
    %       obj.switchDisks = 2;
    %       obj.switchContrast = 0.5; % 50% contrast
    %%%
    
    % Slice disks into pie-shaped regions
    settings.slices              = 1;    % How many pie-shaped regions?
        settings.sliceDisks      = [];   % Which disk #s recieve slices?
        settings.sliceRotation   = 0;    % Where to place disks rotationally? (degrees)
        
    %%% EXAMPLES
    %   - Split center disk into 7 slices
    %       obj.slices = 7;
    %       obj.sliceDisks = 1;
    %   - Split second-most center disk into 3 slices rotated 30 deg.
    %       obj.slices = 3;
    %       obj.sliceDisks = 2;
    %       obj.sliceRotation = 30;
    %%%

    settings.smoothing           = false; % Whether to smooth edges on final movies 
    settings.numberOfMetamerMovies = 1;   % Number of unique metamers to make
    
end