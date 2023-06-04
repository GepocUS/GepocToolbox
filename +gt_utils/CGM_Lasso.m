%% CGM_Lasso - Computes the Composite Gradient Mapping operator for Lasso problem
% 
% This function is used in example_RFISTA
%
% This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function ymas = CGM_Lasso(y, g, R, lambda)

    q = g(y) - R*y;
    s_q = sign(q);

    lambda = s_q.*lambda;
    ymas = R\(-q + lambda);
    pos = find(s_q == 1);
    neg = find(s_q == -1);
    ymas(pos) = min( ymas(pos), zeros(length(pos), 1) );
    ymas(neg) = max( ymas(neg), zeros(length(neg), 1) );

end

