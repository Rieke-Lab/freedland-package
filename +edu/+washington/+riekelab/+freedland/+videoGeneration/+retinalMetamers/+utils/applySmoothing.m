% Function within retinalMetamers
% By J. Freedland, 2020
%
%%% Applies Gaussian smoothing along edges of image. 

function movieOutput = applySmoothing(movieInput,masks)

    width = 5;
    smoothingMask = zeros(size(movieInput,1),size(movieInput,2));
    
    for a = 1:size(masks,3)
        S = masks(:,:,a);

        % Shrink mask in all directions and use padding.
        S_small = imresize(S,[size(S,1)-width,size(S,2)-width],'nearest');
        S_small = [zeros(size(S_small,1),width), S_small, zeros(size(S_small,1),width)];
        S_small = [zeros(width,size(S_small,2)); S_small; zeros(width,size(S_small,2))];

        % Enlarge mask in all directions.
        S_large = imresize(S,[size(S,1)+width,size(S,2)+width],'nearest');

        % Build separated mask
        S2 = abs(S_large - S_small);
        blend = imresize(S2,[size(S,1),size(S,2)]);
        blend = (blend ~=0);
        
        % Add to iterative mask
        smoothingMask = smoothingMask + blend;
    end
    
    smoothingMask(smoothingMask > 1) = 1;

    movieOutput = zeros(size(movieInput));
    for b = 1:size(movieInput,5) % Each movie
        for c = 1:size(movieInput,4) % Each frame

            % Blur relevant frame
            regionBlur = imgaussfilt(movieInput(:,:,:,c,b),2);
            regionBlur = regionBlur .* smoothingMask; % Only take specific region

            % Remove region from movie
            movieOutput(:,:,:,c,b) = movieInput(:,:,:,c,b) .* double(smoothingMask == 0); % Remove region
            
            % Add in blurred region
            movieOutput(:,:,:,c,b) = movieOutput(:,:,:,c,b) + regionBlur; % re-add frame

        end
    end
end

