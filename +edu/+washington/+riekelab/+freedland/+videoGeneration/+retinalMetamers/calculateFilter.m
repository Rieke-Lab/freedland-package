% Function within retinalMetamers
% By J. Freedland, 2020
%
% After calculating a neuron's receptive field (RF), a difference of gaussian
% (DoG) filter is created. This is used to normalize values in the
% naturalistic image's projection.
%%%

function [normFilter,f,f_um] = calculateFilter(obj,normalization)

    % Convert neuron's RF to DOVES VH units.
    centerSigmaPix = retinalMetamers.utils.changeUnits(obj.rfSigmaCenter,obj.micronsPerPixel,'UM2VH');
    surroundSigmaPix = retinalMetamers.utils.changeUnits(obj.rfSigmaSurround,obj.micronsPerPixel,'UM2VH');

    % Generate 2D gaussians
    centerGaus = fspecial('gaussian',[obj.videoSize(1) obj.videoSize(2)],centerSigmaPix);
    surroundGaus = fspecial('gaussian',[obj.videoSize(1) obj.videoSize(2)],surroundSigmaPix);

    % Calculate difference of gaussians
    diffGaussian = centerGaus - surroundGaus;
    normFilter = diffGaussian ./ max(diffGaussian(:)); % Normalize filter

    % Take 2D slice of half-gaussian.
    slice       = normFilter(round(obj.videoSize(1)/2),round(obj.videoSize(2)/2):end); 
    curvature   = diff(slice);      % Derivative
    tot         = cumsum(slice);    
    tot         = tot ./ max(tot(:));% Normalized integral
    
    % Calculate specific regions
    f = cell(7,2);
    
    f{2,1} = 'Maximally inhibitory';
    [~,f{2,2}] = min(slice); % Location with largest inhibitory response.
    
    f{1,1} = 'Zero point';
    [~,f{1,2}] = min(abs(slice(1:f{2,2}))); % Where excitatory response switches to inhibitory.
    
    f{3,1} = 'Largest [excitatory inhibitory] curvature';
    [~,a] = min(curvature);     % Excitatory curvature
    [~,b] = max(curvature);     % Inhibitory curvature
    f{3,2} = [a b];

    f{4,1} = 'Excitatory/inhibitory balance point';
    [~,f{4,2}] = min(abs(tot(1:f{2,2}) - tot(end))); % Where total excitation = total inhibition
    
    f{5,1} = '% Excitatory';
    range = 5:5:95;
    percentage = zeros(1,length(range));
    for a = 1:length(range)
        [~,percentage(a)] = min(abs(tot(1:f{1,2}) - range(a)/100));
    end
    f{5,2} = [range; percentage];
    
    f{6,1} = '% Inhibitory';
    range = 5:5:95;
    percentage = zeros(1,length(range));
    regime = (1 - tot(f{1,2}:end)) / (1-tot(end)); % Only consider inhibitory region
    for a = 1:length(range)
        [~,percentage(a)] = min(abs(regime - range(a)/100));
        percentage(a) = percentage(a) + f{1,2};
    end
    f{6,2} = [range; percentage];
    
    f{7,1} = 'Total % Inhibitory';
    f{7,2} = 1 - tot(end);
    
    % Convert to microns (from pixels)
    f_um = f;
    for a = 1:6
        f_um{a,2} = retinalMetamers.utils.changeUnits(f_um{a,2},obj.micronsPerPixel,'VH2UM');
    end
    f_um{5,2}(1,:) = f{5,2}(1,:); % Percentages do not change
    f_um{6,2}(1,:) = f{6,2}(1,:);

    % Shift gaussian s.t. all values are positive
    if normalization == true
        shift = abs(min(normFilter(:)));  
        normFilter = normFilter + shift;
        normFilter = normFilter ./ max(normFilter(:)); % Renormalize
    end

% %     To visualize this region, uncomment:
%     figure(1)
%     umWidth = retinalMetamers.utils.changeUnits(1:length(slice),obj.micronsPerPixel,'VH2UM');
%     plot(umWidth,slice)
%     hold on
%     plot(f_um{1,2},slice(f{1,2}),'ro','LineWidth',2)
%     plot(f_um{2,2},slice(f{2,2}),'bo','LineWidth',2)
%     plot(f_um{3,2},slice(f{3,2}),'go','LineWidth',2)
%     plot(f_um{4,2},slice(f{4,2}),'yo','LineWidth',2)
%     plot(f_um{5,2}(2,:),slice(f{5,2}(2,:)),'kx')
%     plot(f_um{6,2}(2,:),slice(f{6,2}(2,:)),'kx')
%     legend('RF',f_um{1,1},f_um{2,1},f_um{3,1},f_um{4,1},f_um{5,1},f_um{6,1})
%     hold off
%     xlabel('microns')
%     ylabel('neuron weight')
%     waitforbuttonpress 
end