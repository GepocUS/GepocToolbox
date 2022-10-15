%% PellipProj - Perform P-wighted projection onto an ellipsoid of geometry P
%
% Returns the solution of the optimization problem:
% 
%  v = min_x || x - z ||^2_P
%      s.t.  (x - c)'*P*(x - c) <= r^2
%
% for the given z, P, c and r. It is assumed that r > 0.
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
