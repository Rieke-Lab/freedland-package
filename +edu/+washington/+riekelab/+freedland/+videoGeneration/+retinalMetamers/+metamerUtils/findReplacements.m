% Function within retinalMetamers
% By J. Freedland, 2020
%
% Builds metamers based on a library of similarly projected regions.
%
% OUTPUT:   stimulus.metamer:   collection of metamer movies
%                               a 5D matrix. Each movie spans the first 4
%                               dimensions (i.e.
%                               stimulus.metamer(:,:,:,:,1) is one
%                               generated movie.
%
%                               note: the first movie, stimulus.metamer(:,:,:,:,1),
%                               is always the best fit.
%
%           stimulus.metamerProjection: low-dimensional projection of each
%                                       metamer.
%
%           stimulus.error:     difference (in percent contrast) between
%                               the natural image projection and metamer
%                               projection.
%
%           stimulus.metamerValues: light intensity of each projected disk.
%           
%%%

function stimulus = findReplacements(obj,stimulus,databaseValues,databaseTraj)

    stimulus.metamer = zeros(size(stimulus.projection,1),size(stimulus.projection,2),1,size(stimulus.projection,4),obj.numberOfMetamerMovies);
    stimulus.metamerProjection = zeros(size(stimulus.projection,1),size(stimulus.projection,2),1,size(stimulus.projection,4),obj.numberOfMetamerMovies);
    stimulus.error = zeros(obj.numberOfMetamerMovies,size(stimulus.values,2));
    stimulus.metamerValues = zeros(size(stimulus.values,1),size(stimulus.values,2),obj.numberOfMetamerMovies);
    
    for a = 1:size(stimulus.projection,4)

        % Find best replacement for each disk
        A = stimulus.values(:,a);
        B = abs(A - databaseValues);
        [C,i] = sort(B,2);

        % Remove irrelevant rows
        [~,remv] = unique(i,'rows'); % Ignore duplicated masks
        C = C(sort(remv'),:);
        remv = sort(remv);
        remv(sum(C,2) == 0) = [];    % Remove background areas

        % Build best movies
        for c = 1:obj.numberOfMetamerMovies % Find collections of good frames
            stimulus.error(c,a) = nanmean(C(:,c)) ./ obj.backgroundIntensity * 100; % Calculated as percent contrast
            for b = remv'
                stimulus.metamer(:,:,1,a,c) = stimulus.metamer(:,:,1,a,c) + (databaseTraj(:,:,i(b,c)) .* stimulus.masks(:,:,b));
                stimulus.metamerProjection(:,:,1,a,c) = stimulus.metamerProjection(:,:,1,a,c) + (databaseValues(b,i(b,c)) .* stimulus.masks(:,:,b));
                stimulus.metamerValues(b,a,c) = databaseValues(b,i(b,c));
            end
        end
    end
    
    surroundMask = (sum(stimulus.masks,3) == 0);
    stimulus.metamer = stimulus.metamer + surroundMask .* obj.backgroundIntensity;
    stimulus.metamerProjection = stimulus.metamerProjection + surroundMask .* obj.backgroundIntensity;
end