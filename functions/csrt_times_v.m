%% csrt_times_v - Multiplies the transposed of a matrix in CSR format with a vector
% 
%   y = csrt_times_v(val, col, row_ptr, x) computes the multiplication between the
%   vector x and the transposed of the matrix in the Compressed Sparse Row (CSR) form
%   given by the vectors val, col and row_ptr.
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
% See also full2CSR, csr_times_v
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function y = csrt_times_v(val, col, row_ptr, x)

    % Initialize
    n = length(row_ptr)-1;
    y = zeros(max(col), 1);
    
    % Perform multiplication
    for i = 1:n
        for j = row_ptr(i):row_ptr(i+1)-1
                y(col(j)) = y(col(j)) + val(j)*x(i);
        end
    end

end
