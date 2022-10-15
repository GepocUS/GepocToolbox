%% PellipProj - Perform P-wighted projection onto an ellipsoid of geometry P
%
% Returns the solution of the optimization problem:
% 
%  v = min_x || x - z ||^2_P
%      s.t.  (x - c)'*P*(x - c) <= r^2
%
% for the given z, P, c and r. It is assumed that r > 0.
% 
% INPUTS:
%   - z: Vector to be projected onto the ellipsoid
%   - P: Positive definite matrix defining the geometry of the ellipsoid
%   - c: Vector defining the center of the ellipsoid
%   - r: scalar defining the size of the ellipsoid
%
% OUTPUTS:
%   - v: P-weignted projection onto the ellipsoid of geometry P
%
% This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
%
    
function v = PellipProj(z, P, c, r)

    % Compute (z - c)'*P*(z - c)
    zPz = (z - c)'*P*(z - c);
    
    % Compute projection
    if zPz <= r^2
        v = z; % z belongs to the ellipsoid
    else
        v = r*(z - c)/sqrt(zPz) + c;
    end

end
