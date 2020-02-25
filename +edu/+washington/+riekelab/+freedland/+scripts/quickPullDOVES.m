function [path, outputPicture] = quickPullDOVES(imageNo, observerNo, mirroring)
    
    % Format imageNumber into correct filename.
    imageString = num2str(imageNo);
    directory = '+edu/+washington/+riekelab/+freedland/+images/';
    
    while size(imageString,2) < 3
        imageString = strcat('0',imageString);
    end
    
    filename = strcat(directory,'img',imageString,'.mat');
    
    load(filename)
    
    % Assign data to easy-to-use variables.
    observer = eye_data{1,observerNo};
    path.x = round(observer(1,:)); % Assign x coordinate data
    path.y = round(observer(2,:)); % Assign y coordinate data
    
    % Flip if needed
    if mirroring == true
        mirroredImage = [flip(flip(picture,2)) flip(picture) flip(flip(picture,2));
        flip(picture,2) picture flip(picture,2);
        flip(flip(picture,2)) flip(picture) flip(flip(picture,2))];
    
        outputPicture = mirroredImage;
    
        c = size(picture);
        
        path.x = path.x + c(2);
        path.y = path.y + c(1);
    else
        outputPicture = picture;
    end
    
end