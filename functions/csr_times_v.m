%% csr_times_v - Multiplies matrix in CSR format with a vector
% 
%   y = csr_times_v(val, col, row_ptr, x) computes the multiplication between the
%   vector x and the matrix in the Compressed Sparse Row (CSR) form given
%   by the vectors val, col and row_ptr.
%
% INPUTS:
%   - var: Vector containing the values of the non-zero elements.
%   - col: Contains the column-index of each of the non-zero elements in var.
%   - row_ptr: Contains the row-index range of the non-zero elements in var.
%   - x: Vector to which to multiply the matrix by.
%
% OUTPUS:
%   - y: Result of the matrix-vector multiplication.
%
% See also full2CSR.
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function y = csr_times_v(val, col, row_ptr, x)

    % Initialize
    n = length(row_ptr)-1;
    y = zeros(n, 1);
    
    % Perform multiplication
    for i = 1:n
        for j = row_ptr(i):row_ptr(i+1)-1
                y(i) = y(i) + val(j)*x(col(j));
        end
    end

end
