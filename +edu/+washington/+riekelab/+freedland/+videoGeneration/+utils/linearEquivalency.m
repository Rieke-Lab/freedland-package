% Function within retinalMetamers
% By J. Freedland, 2020
%
% The head honcho. Calculates weighted averages (using a natural image
% movie convolved with the neuron's RF) for each pre-specifed region.
%
% OUTPUT:   newTraj: projected naturalistic image.
%           diskValues: light intensity of each disk
%           masks: shape of each mask
%%%
function [newTraj, diskValues, masks] = linearEquivalency(retinalMetamers, weightedTrajectory, RFFilter, unweightedTrajectory)   

    % Define space as polar coordinates (r = radial, th = theta)
    [xx,yy] = meshgrid(1:retinalMetamers.videoSize(2),1:retinalMetamers.videoSize(1));
    r = sqrt((xx - retinalMetamers.videoSize(2)/2).^2 + (yy - retinalMetamers.videoSize(1)/2).^2); 
    th = atan((xx - retinalMetamers.videoSize(2)/2) ./ (yy - retinalMetamers.videoSize(1)/2));
    th = abs(th-pi/2);              
    
    % Adjust theta space for strange monitors
    nonsmooth = find(diff(th) > pi/2,1);
    th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
    th = rad2deg(th);
    th = mod(th + retinalMetamers.sliceRotation,360); % Rotate as required
    
    % Identify user-specific radii
    if strcmp(retinalMetamers.diskRegionUnits,'pix')
        retinalMetamers.radii = round(utils.changeUnits(retinalMetamers.diskRegions,retinalMetamers.micronsPerPixel,'pix2arcmin'));
    elseif strcmp(retinalMetamers.diskRegionUnits,'um')
        retinalMetamers.radii = round(utils.changeUnits(retinalMetamers.diskRegions,retinalMetamers.micronsPerPixel,'um2arcmin'));
    elseif strcmp(retinalMetamers.diskRegionUnits,'deg')
        retinalMetamers.radii = retinalMetamers.radii / 60; % DOVES units are in arcmin
    elseif strcmp(retinalMetamers.diskRegionUnits,'arcmin')
        retinalMetamers.radii = retinalMetamers.diskRegions;
    else
        error('Please identify correct diskRegionUnit: "pix", "um", "deg", or "arcmin"')
    end
    
    retinalMetamers.slices(retinalMetamers.slices == 0) = 1;
    retinalMetamers.theta = 0:360/retinalMetamers.slices:360;
    
    % Generate switchDisk values, used to measure impact of surround.
    if sum(retinalMetamers.switchDisks) > 0
        retinalMetamers.switchTraj = utils.generateSwitchDisks(retinalMetamers);
    end

    newTraj = zeros(size(weightedTrajectory));
    diskValues = zeros((length(retinalMetamers.radii)-1).*(length(retinalMetamers.theta)-1),size(weightedTrajectory,4));
    masks = zeros(size(weightedTrajectory,1),size(weightedTrajectory,2),(length(retinalMetamers.radii)-1).*(length(retinalMetamers.theta)-1));

    % Calculate statistics
    counter = 1;
    for a = 1:length(retinalMetamers.radii) - 1
        radiusFilt = r >= retinalMetamers.radii(a) & r <= retinalMetamers.radii(a+1); % Radial filter (r)
        for b = 1:length(retinalMetamers.theta) - 1
            angFilt = th >= retinalMetamers.theta(b) & th < retinalMetamers.theta(b+1); % Angular filter (theta)
            ignoreDisk = false;
            
            % Whether to ignore angular filter
            if ~ismember(a,retinalMetamers.sliceDisks)
                angFilt = ones(size(angFilt));
                if b > 1
                    ignoreDisk = true;
                end
            end  

            filt = radiusFilt .* angFilt;
            tempMask = zeros(size(weightedTrajectory));
            
            % For linear equivalent regions
            if ignoreDisk == false
                masks(:,:,counter) = filt;
                if ismember(a,retinalMetamers.meanDisks) || ismember(a,retinalMetamers.metamerDisks)
                    
                    % Normalizing value
                    T = RFFilter .* filt;
                    T = sum(T(:));

                    % Apply across trajectory
                    for c = 1:size(weightedTrajectory,4)
                        S = weightedTrajectory(:,:,1,c) .* filt;
                        diskValues(counter,c) = sum(S(:)) / T;
                        tempMask(:,:,1,c) = diskValues(counter,c) .* filt;
                    end
                elseif ismember(a,retinalMetamers.backgroundDisks)
                    
                    % Apply across trajectory
                    for c = 1:size(diskValues,2)
                        tempMask(:,:,1,c) = retinalMetamers.backgroundIntensity .* filt;
                        diskValues(counter,c,:) = retinalMetamers.backgroundIntensity;
                    end
                elseif ismember(a,retinalMetamers.naturalDisks)
                    
                    % Apply original image to region
                    for c = 1:size(diskValues,2)
                        tempMask(:,:,1,c) = unweightedTrajectory(:,:,1,c) .* filt;
                        diskValues(counter,c,:) = NaN;
                    end

                elseif ismember(a,retinalMetamers.switchDisks)
                    % Apply specific intensity to region
                    for c = 1:size(diskValues,2)
                        tempMask(:,:,1,c) = retinalMetamers.switchTraj(c) .* filt;
                        diskValues(counter,c,:) = retinalMetamers.switchTraj(c);
                    end
                    
                else
                    error(strcat('please specify disk type for disk #',num2str(a)))
                end
                newTraj = newTraj + tempMask;
            end
            counter = counter + 1;
        end
    end

    % Add surrounding mask
    surroundMask = (sum(masks,3) == 0);
    newTraj = newTraj + repmat(surroundMask .* retinalMetamers.backgroundIntensity,1,1,1,size(newTraj,4));

    % Remove excess values
    rm = find(squeeze(sum(sum(masks,1),2)) == 0);
    diskValues(rm,:,:) = [];
    masks(:,:,rm) = [];
end