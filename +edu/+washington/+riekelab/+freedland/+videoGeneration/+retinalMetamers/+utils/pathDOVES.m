% Function within retinalMetamers
% By J. Freedland, 2020
%
%%% Pulls information from DOVES database.

function [path, im, saccades] = pathDOVES(imageNo, observerNo) 

    % Load DOVES information from directory
    imageString = num2str(imageNo);
    directory = '+edu/+washington/+riekelab/+freedland/+images/';
    while size(imageString,2) < 3
        imageString = strcat('0',imageString);
    end
    load(strcat(directory,'img',imageString,'.mat'))

    % Pull relevant information
    observer = eye_data{1,observerNo};
    path.x = observer(1,:); % x coordinate data
    path.y = observer(2,:); % y coordinate data
        
    % Identify saccades
    movement = diff(path.x).^2 + diff(path.y).^2; % total eye movement
    saccades = find(movement > 5) + 1;  % 60px = 1 deg

    % Mirror image to prevent clipping
    im = [flip(flip(picture,2)) flip(picture) flip(flip(picture,2));
    flip(picture,2) picture flip(picture,2);
    flip(flip(picture,2)) flip(picture) flip(flip(picture,2))];
    path.x = round(path.x + size(picture,2));
    path.y = round(path.y + size(picture,1));
end