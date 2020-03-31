% Function within retinalMetamers
% By J. Freedland, 2020
%
% After calculating a neuron's receptive field (RF), a difference of gaussian
% (DoG) filter is created. This is used to normalize values in the
% naturalistic image's projection.
%%%

function [normFilter,f] = calculateFilter(obj)

    % Convert neuron's RF to DOVES VH units.
    centerSigmaPix = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.rfSigmaCenter,obj.micronsPerPixel,'UM2VH');
    surroundSigmaPix = edu.washington.riekelab.freedland.videoGeneration.retinalMetamers.utils.changeUnits(obj.rfSigmaSurround,obj.micronsPerPixel,'UM2VH');

    % Generate 2D gaussians
    centerGaus = fspecial('gaussian',[obj.videoSize(1) obj.videoSize(2)],centerSigmaPix);
    surroundGaus = fspecial('gaussian',[obj.videoSize(1) obj.videoSize(2)],surroundSigmaPix);

    % Calculate difference of gaussians
    diffGaussian = centerGaus - surroundGaus;
    normFilter = diffGaussian ./ max(diffGaussian(:)); % Normalize filter

    % Shift gaussian s.t. all values are positive
    shift = abs(min(normFilter(:)));  
    normFilter = normFilter + shift;
    normFilter = normFilter ./ max(normFilter(:)); % Renormalize

    % Take 2D slice of half-gaussian.
    slice       = normFilter(round(obj.videoSize(1)/2),round(obj.videoSize(2)/2):end);
    curvature   = diff(slice);
    
    % Use inhibitory region to find symmetrical regions
    f = zeros(2,1);
    inh = abs(find(curvature > 0, 1));
    [~,f(2)] = max(curvature);
    [~,f(1)] = min(abs(slice(1:inh) - slice(f(2,1))));

    % To visualize this region, uncomment:
%     figure(1)
%     plot(slice)
%     hold on
%     plot(f(1):f(2),slice(f(1):f(2)),'r')
%     hold off
%     xlabel('# of pixels from center')
%     ylabel('neuron weight')
%     waitforbuttonpress
end