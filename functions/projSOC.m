%% projSOC - projection to the second order cone (SOC)
%
% Returns the projection of z \in \R^n onto the SOC
%
% SOC = { x = (x_0, x_1) \in \R \times \R^{n-1} : || x_1 ||_2 <= x_0 }
%
% INPUTS:
%   - z: Point to project onto the SOC
%
% OUTPUTS:
%   - v: Projection of z onto the SOC
%
% This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
%

function v = projSOC(z)

    %% Extract z_0, z_1 and compute norm of z_1
    n = length(z);
    z_0 = z(1);
    z_1 = z(2:end);
    nz_1 = norm(z_1, 2);
    
    %% Project
    if nz_1 <= z_0
        v = z;
        
    elseif nz_1 <= -z_0
        v = zeros(n, 1);
        
    else
        v = ((z_0 + nz_1)/(2*nz_1))*[nz_1; z_1];
    end

end
