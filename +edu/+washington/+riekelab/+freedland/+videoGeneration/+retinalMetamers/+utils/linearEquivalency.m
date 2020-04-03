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
function [newTraj, diskValues, masks] = linearEquivalency(obj, weightedTrajectory, RFFilter, unweightedTrajectory)   
    
    % Define space as polar coordinates (r = radial, th = theta)
    [xx,yy] = meshgrid(1:obj.videoSize(2),1:obj.videoSize(1));
    r = sqrt((xx - obj.videoSize(2)/2).^2 + (yy - obj.videoSize(1)/2).^2); 
    th = atan((xx - obj.videoSize(2)/2) ./ (yy - obj.videoSize(1)/2));
    th = abs(th-pi/2);              
    
    % Adjust theta space for strange monitors
    nonsmooth = find(diff(th) > pi/2,1);
    th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
    th = rad2deg(th);
    
    % Identify user-specific radii
    if strcmp(obj.diskRegionUnits,'pixels')
        obj.radii = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.diskRegions,obj.micronsPerPixel,'PIX2VH'));
    elseif strcmp(obj.diskRegionUnits,'microns')
        obj.radii = round(edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.diskRegions,obj.micronsPerPixel,'UM2VH'));
    elseif strcmp(obj.diskRegionUnits,'deg')
        obj.radii = obj.radii / 60; % DOVES units are in arcmin
    elseif strcmp(obj.diskRegionUnits,'rf')
        obj.radii = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.RFConversion(obj);
    elseif strcmp(obj.diskRegionUnits,'arcmin')
        obj.radii = obj.diskRegions;
    else
        error('incorrect diskRegionUnit')
    end
    
    if obj.slices > 0
        obj.theta = 0:360/obj.slices:360;
    else
        obj.theta = [0 360];
    end
    
    % Generate switchDisk values, used to measure impact of surround.
    if sum(obj.switchDisks) > 0
        obj.switchTraj = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.generateSwitchDisks(obj);
    end

    diskValues = zeros(length(obj.radii-1).*length(obj.theta)-1,size(weightedTrajectory,4));
    newTraj = zeros(size(weightedTrajectory));
    masks = zeros(size(weightedTrajectory,1),size(weightedTrajectory,2),size(diskValues,1));

    % Calculate statistics
    counter = 1;
    for a = 1:length(obj.radii) - 1
        radiusFilt = r >= obj.radii(a) & r <= obj.radii(a+1); % Radial filter (r)
        for b = 1:length(obj.theta) - 1
            angFilt = th >= obj.theta(b) & th < obj.theta(b+1); % Angular filter (theta)
            
            if ~ismember(a,obj.cuts) % Ignore angular filter
                angFilt = ones(size(angFilt));
            end  

            filt = radiusFilt .* angFilt;
            masks(:,:,counter) = filt;

            tempMask = zeros(size(weightedTrajectory));
            
            % For linear equivalent regions
            if ismember(a,obj.meanDisks)
                
                % Normalizing value
                regionSize = sum(filt(:));
                T = RFFilter .* filt;
                R = sum(T(:)) / regionSize;

                % Apply across trajectory
                for c = 1:size(weightedTrajectory,4)
                    S = weightedTrajectory(:,:,1,c) .* filt;
                    
                    diskValues(counter,c) = sum(S(:)) / (regionSize * R);
                    tempMask(:,:,1,c) = diskValues(counter,c) .* filt;
                end

            elseif ismember(a,obj.backgroundDisks)

                % Apply across trajectory
                for c = 1:size(diskValues,2)
                    tempMask(:,:,1,c) = obj.backgroundIntensity .* filt;
                    diskValues(counter,c) = obj.backgroundIntensity;
                end

            elseif ismember(a,obj.switchDisks)

                % Apply specific intensity to region
                for c = 1:size(diskValues,2)
                    tempMask(:,:,1,c) = obj.switchTraj(c) .* filt;
                    diskValues(counter,c) = obj.switchTraj(c);
                end
                
            elseif ismember(a,obj.naturalDisks)
                
                % Apply original image to region
                for c = 1:size(diskValues,2)
                    tempMask(:,:,1,c) = unweightedTrajectory(:,:,1,c) .* filt;
                    diskValues(counter,c) = NaN;
                end

            end

            if ~(~ismember(a,obj.cuts) && b > 1)
                newTraj = newTraj + tempMask;
            end

            counter = counter + 1;

        end
    end
    
    % Add surrounding mask
    surroundMask = (sum(masks,3) == 0);
    newTraj = newTraj + surroundMask .* obj.backgroundIntensity;
end