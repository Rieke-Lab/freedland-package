% We use the DOVES database to produce an incident image with and without
% fixational eye movements, or drift.
%
% INPUT:    imageNo: natural image from DOVES database, #1-101.
%           observerNo: specific observer from DOVES database. #1-19.
%   varargin:       'offSetWidth', #: offset from focused point (in degrees)
%                   'offSetHeight', #: offset from focused point (in degrees)
%                   'show...', true: plays resulting movie.
%
% OUTPUT:   empiricalMotion: video with both saccades and drift.
%           saccadeMotion: video with only saccades. (drift removed).
%           randomMotion: video with saccades + random walk.
%           information: contains original image and helpful information for analysis.
%
% Created by JMF 09/2018
% Last edited by JMF 09/2018

function [empiricalMotion, saccadeMotion, randomMotion, information] = arrangeDOVES(imageNo, observerNo, varargin) 
    
    ip = inputParser(); % Used to dictate settings.
    ip.addRequired('imageNo', @isnumeric);
    ip.addRequired('observerNo', @isnumeric);
    
    ip.addParameter('sigma', 3, @isnumeric); % Standard deviation of random walk in microns
    ip.addParameter('directory', true, @islogical);
    % If images are in the same folder, make false. If not, type path in 'directory' variable.
    
    ip.addParameter('offSetHeight', 0, @isnumeric); % Offset from center of eye in microns
    ip.addParameter('offSetWidth', 0, @isnumeric);
    
    ip.addParameter('pixelHeight', 17, @isnumeric); % Pixels in field of view
    ip.addParameter('pixelWidth', 17, @isnumeric);
    
    ip.addParameter('adjustResolution', false, @islogical);
    ip.addParameter('linearlyAdjust', [600 1600], @isnumeric);
    
    ip.parse(imageNo, observerNo, varargin{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Format imageNumber into correct filename.
    imageString = num2str(imageNo);
    directory = '+edu/+washington/+riekelab/+freedland/+images/';
    
    while size(imageString,2) < 3
        imageString = strcat('0',imageString);
    end
    if ip.Results.directory
        filename = strcat(directory,'img',imageString,'.mat');
    else
        filename = strcat('img',imageString,'.mat');
    end
    
    load(filename)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign variables
    imHeight = ip.Results.pixelHeight-1;
    imWidth = ip.Results.pixelWidth-1;
    
    if mod(imHeight,2) == 1 % Make even values odd
        imHeight = imHeight - 1;
    end
    if mod(imWidth,2) == 1 % Make even values odd
        imWidth = imWidth - 1;
    end
    
    % Assign data to easy-to-use variables.
    observer = eye_data{1,observerNo};
    observerX = observer(1,:); % Assign x coordinate data
    observerY = observer(2,:); % Assign y coordinate data
    offSetHeight = ip.Results.offSetHeight; % Assign offsets
    offSetWidth = ip.Results.offSetWidth; % Assign offsets
    reach = [floor(imWidth/2) floor(imHeight/2)]; % Defines range to analyze
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % We want our retina to always see an image, even at corners. This will
    % produce a mirrored effect if we reach an edge.
    mirroredImage = [flip(flip(picture,2)) flip(picture) flip(flip(picture,2));
    flip(picture,2) picture flip(picture,2);
    flip(flip(picture,2)) flip(picture) flip(flip(picture,2))];

    % And create a coordinate offset to look around the mirrored image
    c = size(picture);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Pre-allocate memory.
    xDev = zeros(1,size(observer,2));
    yDev = zeros(1,size(observer,2));
    empirical = zeros(imHeight+1,imWidth+1,1,size(observer,2));
    
    % Build n x m x 1 x a matrix where 'a' contains the image (nxm) at every moment.
    % The third dimension (of size 1) tells MATLAB it is a greyscale image.
    for a = 1:size(observer,2)
        
        imCenter = [round(observerX(1,a))+offSetWidth round(observerY(1,a))+offSetHeight]; % [x y]
        imX = imCenter(1,1) - reach(1,1) : imCenter(1,1) + reach(1,1);
        imY = imCenter(1,2) - reach(1,2) : imCenter(1,2) + reach(1,2);

        % Arrange into matrix
        empirical(:,:,1,a) = mirroredImage(imY+c(1),imX+c(2));
        
        % Important for removing drift in the next stage.        
        xDev(1,a) = median(imX+c(2));
        yDev(1,a) = median(imY+c(1));
    end
    
    saccade = empirical; % Make the matrix. We will then remove drift.
    random = empirical; % This will introduce random new drift.
   
    xChange = zeros(1,size(observer,2)-1);
    yChange = zeros(1,size(observer,2)-1);
    saccadeTracking = zeros(1,size(observer,2));
    
    counter = 1; % Keeps track of when a saccade occurs.
    
    for a = 1:size(observer,2)-1

        % Determines when pixels (in x and y directions) move drastically.
        xChange(1,a) = abs(xDev(1,a) - xDev(1,a+1));
        yChange(1,a) = abs(yDev(1,a) - yDev(1,a+1));
        
        saccadeCutoff = 1/2; % distinguishes how many degrees define a saccade

        % If large movement, replace prior images with a still image.
        if sqrt((xChange(1,a)^2 + yChange(1,a))^2) > saccadeCutoff*60 % in original paper, 60px = 1 deg
            
            saccadeTracking(1,a+1) = 1; % Tracks when a saccade takes place, for later analysis.
            frames = a-counter+1; % number of frames to keep still
            
            xx = round(median(xDev(1,counter:a))) - reach(1,1) : round(median(xDev(1,counter:a))) + reach(1,1); % Use median image
            yy = round(median(yDev(1,counter:a))) - reach(1,2) : round(median(yDev(1,counter:a))) + reach(1,2);
            
            saccade(:,:,1,counter:a) = repelem(mirroredImage(yy,xx),1,1,1,frames); % Overwrite fixational motions.
            
            % Define our random motion
            sigma = obj.rig.getDevice('Stage').um2pix(ip.Results.sigma); % Convert to pixels
            
            deltaX = round(normrnd(0, sigma, frames, 1));
            deltaY = round(normrnd(0, sigma, frames, 1));
            
            XX = xx+deltaX;
            YY = yy+deltaY;
            
            counter2 = 1; % another counter
            
            for q = counter:a
                random(:,:,1,q) = mirroredImage(YY(counter2,:),XX(counter2,:));
                counter2 = counter2 + 1;
            end
            
            counter = a;
        end
    end
    
    % Fill in the final frames of the movie we missed in our loop.
    missedFrames = size(saccade,4)-counter;
    
    D = round(missedFrames / 2) + counter; 
    xx = round(xDev(1,D)) - reach(1,1) : round(xDev(1,D)) + reach(1,1);
    yy = round(yDev(1,D)) - reach(1,2) : round(yDev(1,D)) + reach(1,2);
    saccade(:,:,counter:end) = repelem(mirroredImage(yy,xx),1,1,missedFrames+1);
    
    deltaX = round(normrnd(0, sigma, missedFrames+1, 1));
    deltaY = round(normrnd(0, sigma, missedFrames+1, 1));
    
    XX = xx+deltaX;
    YY = yy+deltaY;
    
    counter2 = 1;
    
    for q = counter:size(random,4)
         random(:,:,1,q) = mirroredImage(YY(counter2,:),XX(counter2,:));
         counter2 = counter2 + 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adjust resolution to desired monitor size
    if ip.Results.adjustResolution
        adjResolution = ip.Results.linearlyAdjust;
        
        % Resize images - this is a slow step
        empiricalMotion.img = imresize(uint8(empirical),adjResolution);
        saccadeMotion.img = imresize(uint8(saccade),adjResolution);
        randomMotion.img = imresize(uint8(random),adjResolution);
        information.picture = uint8(picture);
        information.micronsPerPixel = [adjResolution(1)/imHeight adjResolution(2)/imWidth];
    else 
        empiricalMotion.img = uint8(empirical);
        saccadeMotion.img = uint8(saccade);
        randomMotion.img = uint8(random);
        information.picture = uint8(picture); % save image
    end

    cmap = colormap(gray); % Assign colormap as B&W

    empiricalMotion.mov = immovie(empiricalMotion.img,cmap);
    saccadeMotion.mov = immovie(saccadeMotion.img,cmap);
    randomMotion.mov = immovie(randomMotion.img,cmap);
    
    close
    
    % Fill in with our information
    information.saccadeTracking = saccadeTracking; % Tracks when a saccade is made
    information.randomWalkStD_inMicrons = ip.Results.sigma; % Chosen parameter
    information.offSetHeight_inMicrons = ip.Results.offSetHeight; % Chosen parameter
    information.offSetWidth_inMicrons = ip.Results.offSetWidth; % Chosen parameter
    information.imageName = filename;
    information.observer = observerNo;
    
end