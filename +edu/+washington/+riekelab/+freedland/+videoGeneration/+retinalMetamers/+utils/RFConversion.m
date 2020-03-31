% Function within retinalMetamers
% By J. Freedland, 2020
%
%%% RF coordinates are a special form of units in this project.
%%% They define radial regions (for disks) invariant of a neuron's RF.
%
% Coordinates are defined as follows:
%       0: center of the cell's RF
%       1: start of strongly inhibitory region.
%       2: end of strongly inhibitory region.
%       3: edge of monitor
%
% Regions between these regions are linearly interpolated. i.e:
%       0.75: falls 75% between the center (0) and the start of the 
%             strongly inhibitory region (1).
%       1.50: falls halfway between the start of the strongly inhibitory (1)
%             region and end of the strongly inhibitory region (2)
%       1.00: falls directly on the start of the strongly inhibitory (1)
%       region.
%
% This coordinate system is useful because it is invariant of neuron size.
%%%

function radii = RFConversion(obj)

    centerSize = obj.rfSizing(1);
    annulusSize = obj.rfSizing(2) - obj.rfSizing(1);
    surroundSize = (max(obj.videoSize) / 2) - obj.rfSizing(2);

    radii = zeros(obj.numberOfDisks+1,1);
    for a = 1:obj.numberOfDisks+1
        if obj.diskRegions(a) <= 1
            RFCoordinate = obj.diskRegions(a);
            radii(a) = RFCoordinate .* centerSize;
        elseif obj.diskRegions(a) <= 2
            RFCoordinate = obj.diskRegions(a) - 1;
            radii(a) = centerSize + (RFCoordinate .* annulusSize);
        elseif obj.diskRegions(a) <= 3
            RFCoordinate = obj.diskRegions(a) - 2;
            radii(a) = centerSize + annulusSize + (RFCoordinate .* surroundSize);
        end
    end
end