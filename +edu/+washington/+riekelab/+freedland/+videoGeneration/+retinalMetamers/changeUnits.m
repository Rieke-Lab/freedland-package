% Function within retinalMetamers
% By J. Freedland, 2020
%
%%% Change units within retinalMetamers.
%
% Using user settings for "micronsPerPixel" (um/pix), we define a series of units:
% VH:   units used in eye tracking studies, courtesy of the DOVES database
%       (200Hz sampling)
% UM:   microns (across the retina)
% PIX:  monitor pixels (output).
%
% micronsPerPixel:  number of microns each output monitor pixel spans.
%                   depends on experimental setup.
%
% Created by JMF 03/2020
% Last edited by JMF 03/2020
%%%

function B = changeUnits(A,micronsPerPixel,type)
            
    if strcmp(type,'UM2PIX')
        % um / (um/pix) = pix
        B = A ./ micronsPerPixel;

    elseif strcmp(type,'PIX2UM')
        % pix * (um/pix) = um
        B = A .* micronsPerPixel;

    elseif strcmp(type,'UM2VH')
        % From DOVES database: 1 VH pixel = 1 arcmin = 3.3 um on monkey retina
        % um / (3.3 um/VH) = VH
        B = A ./ 3.3;

    elseif strcmp(type,'VH2UM')
        % VH * (3.3 um/VH) = um
        B = A .* 3.3;

    elseif strcmp(type,'PIX2VH')
        % (3.3 um/VH) / (um/pix) = pix/VH
        ratio = 3.3 ./ micronsPerPixel;

        % pix / (pix/VH) = VH
        B = A ./ ratio;

    elseif strcmp(type,'VH2PIX')
        % (3.3 um/VH) / (um/pix) = pix/VH
        ratio = 3.3 ./ micronsPerPixel;

        % VH * (pix/VH) = pix
        B = A .* ratio;
    else
        error('incorrect unit conversion.')
    end
end