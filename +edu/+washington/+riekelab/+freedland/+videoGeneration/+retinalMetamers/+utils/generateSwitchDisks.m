% Function within retinalMetamers
% By J. Freedland, 2020
%
% Generates a "switchDisk" trajectory. These types of disks use saccades to
% define when to switch between bright and dark regions (relative to the
% average light intensity).
%%%
function traj = generateSwitchDisks(obj)

    % Downsample trajectory
    obj.saccades = obj.saccades / (200 / obj.monitorFrameRate);

    % Identify frames to switch disk intensity
    switchLocations = (diff(obj.saccades) > (2 * obj.monitorFrameRate / 200)) .* obj.saccades(2:end);
    switchLocations(switchLocations == 0) = [];
    switchLocations = switchLocations + 1; % Occur as saccade occurs
    switchLocations = round([1 switchLocations length(obj.xTraj)]);

    % Identify user preference
    if strcmp(obj.switchPref,'darkFirst')
        intensities = [0 obj.backgroundIntensity.*2]; % [min, max] intensity between disks
    elseif strcmp(obj.switchPref,'brightFirst')
        intensities = [obj.backgroundIntensity.*2 0]; % swap order
    else
        error('Please specify either darkFirst or brightFirst for obj.switchPref')
    end

    % Build trajectory
    traj = zeros(1,length(obj.xTraj));
    counter = 0;
    for a = 1:length(switchLocations)-1
        traj(switchLocations(a):switchLocations(a+1)) = intensities(counter+1);
        counter = mod(counter+1,2);
    end
end