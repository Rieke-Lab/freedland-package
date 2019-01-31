% This function uses the output of pathDOVES to calculate local image
% statistics across a trajectory.
% INPUT:    aperature diameter: center diameter in px
%           surr diameter: surrounding diameter in px
%           integration function: 'gaussian' or other.
% OUTPUT:   structures with various 1xn vectors of local statistics over a
% trajectory.
% Created: JMF 11.2018

function [center, annulus, surround] = imageSpotMean(xTraj,yTraj,img,canvasSize,RFFilter,obj,varargin)
    
    ip = inputParser(); % Used to dictate settings.
    ip.addRequired('xTraj');
    ip.addRequired('yTraj');
    ip.addRequired('img');
    ip.addRequired('canvasSize')
    ip.addRequired('RFFilter');
    ip.addRequired('obj');
    ip.addParameter('linearEquiv', true);
    ip.addParameter('umPerPix',1.2);
    ip.parse(xTraj, yTraj, img, canvasSize, RFFilter, obj, varargin{:});

    % create distance matrix
    [xx, yy] = meshgrid(1:2*canvasSize(2),1:2*canvasSize(1));
    m = sqrt(((xx-(canvasSize(2))).*obj.compressHorizontally).^2+...
        ((yy-(canvasSize(1))).*obj.compressVertically).^2);
    
    % make elliptical
    r = imrotate(m,obj.rotation);
    rm = r(round(size(r,1)/2 - canvasSize(2)/2):round(size(r,1)/2 + canvasSize(2)/2),...
                round(size(r,2)/2 - canvasSize(1)/2):round(size(r,2)/2 + canvasSize(1)/2));
    
    center.filt = rm < obj.centerRadius.maskR * (1+(obj.expandCenter/100));
    annulus.filt = rm > obj.centerRadius.maskR & rm < obj.surroundRadius.maskR;
    surround.filt = rm > (obj.surroundRadius.maskR * (1-(obj.contractSurround/100)));
    
    center.gradient = double(center.filt) .* RFFilter;
    annulus.gradient = double(annulus.filt) .* RFFilter;
    surround.gradient = double(surround.filt) .* RFFilter;
    
    % we will adjust the surround so it extends beyond our screen.
    % this prevents pixels from "leaking"
    
    surround.safetyVal = [ones(size(surround.filt,1),canvasSize(1)/2) surround.filt...
        ones(size(surround.filt,1),canvasSize(1)/2)];
    surround.safety = [ones(canvasSize(2)/2,size(surround.safetyVal,2));surround.safetyVal;...
        ones(canvasSize(2)/2,size(surround.safetyVal,2))];
    
    % and do the same for our gradient:
    val = surround.gradient(canvasSize(2),canvasSize(1)); % last value
    surround.safetyGrad = [ones(size(surround.filt,1),canvasSize(1)/2).*val surround.gradient...
        ones(size(surround.filt,1),canvasSize(1)/2).*val];
    surround.safetyGradient = [ones(canvasSize(2)/2,size(surround.safetyGrad,2)).*val;surround.safetyGrad;...
        ones(canvasSize(2)/2,size(surround.safetyGrad,2)).*val];
    
    if ip.Results.linearEquiv == true
        % this will appear fairly complicated: in brief, we have sized out
        % our filters successfully. However, we will need to predict
        % our entire trajectory in VH (monitor units) to calculate our
        % linear equivalent disc.
        
        % Based on the DOVES trajectory, we introduce a special parameter.
        VHConv = 3.3/ip.Results.umPerPix;
        
        % all of our parameters remain the same, except for the size
        linEqCanvasSize = ceil(canvasSize ./ VHConv);
        linEqSpace = rm((canvasSize(2)/2+1)-floor(linEqCanvasSize(2)/2):(canvasSize(2)/2+1)+floor(linEqCanvasSize(2)/2),...
            (canvasSize(1)/2+1)-floor(linEqCanvasSize(1)/2):(canvasSize(1)/2+1)+floor(linEqCanvasSize(1)/2));
        linEqSpace = imresize(linEqSpace, [linEqCanvasSize(2),linEqCanvasSize(1)]); % ensure similar output
        linEqRF = RFFilter((canvasSize(2)/2+1)-floor(linEqCanvasSize(2)/2):(canvasSize(2)/2+1)+floor(linEqCanvasSize(2)/2),...
            (canvasSize(1)/2+1)-floor(linEqCanvasSize(1)/2):(canvasSize(1)/2+1)+floor(linEqCanvasSize(1)/2));
        linEqRF = imresize(linEqRF, [linEqCanvasSize(2),linEqCanvasSize(1)]);
        
        % resize our calculations
        linEqSpaceC = linEqSpace < (obj.centerRadius.linR ./ VHConv);
        linEqFilterC = double(linEqSpaceC) .* linEqRF;
        
        linEqSpaceA = linEqSpace > (obj.centerRadius.maskR./ VHConv) ...
            & linEqSpace < (obj.surroundRadius.maskR./ VHConv);
        linEqFilterA = double(linEqSpaceA) .* linEqRF;
        
        linEqSpaceS = linEqSpace > (obj.surroundRadius.maskR./ VHConv);
        linEqFilterS = double(linEqSpaceS) .* linEqRF;

        for a = 1:size(xTraj,2)
            temp = img(yTraj(1,a)-floor(linEqCanvasSize(2)/2):yTraj(1,a)+floor(linEqCanvasSize(2)/2),...
                xTraj(1,a)-floor(linEqCanvasSize(1)/2):xTraj(1,a)+floor(linEqCanvasSize(1)/2));  

            centerStat = temp .* linEqFilterC;
            centerStat(centerStat == 0) = [];

            surroundStat = temp .* linEqFilterS;
            surroundStat(surroundStat == 0) = [];
            
            annulusStat = temp .* linEqFilterA;
            annulusStat(annulusStat == 0) = [];

            center.data(1,a) = mean(centerStat);
            surround.data(1,a) = mean(surroundStat);
            annulus.data(1,a) = mean(annulusStat);
        end
    end
    
end