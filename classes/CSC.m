%% CSC - Class for representing and working with sparse matrices
% in the Compressed Sparse Column format
%
% TODO:
%   - Add to_dense() method
%   - Add operator overloads (multiplication, addition, transposed, etc.). The same ones that are available for CSR
% 
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox

classdef CSC < GT_SparseMatrix

    methods

    %% CONSTRUCTOR

    function self = CSC(M, threshold)
        % CSR constructor
        % INPUTS:
        %   - M: Matrix to be converted to CSC format
        %   - threshold: Any number smaller than this value is considered 0.
        %                This in an optional argument that defaults to 0.
        if nargin == 0
            return;
        elseif nargin == 1
            threshold = 0; % Any number smaller than threshold is considered a 0
        end

        [~, self.val, self.row, self.col, self.nnz, self.nrow, self.ncol] = CSC.from_dense(M, threshold);

    end

    end % Methods

    methods (Static)

    function [M_sparse, val, row, col, nnz, n, m] = from_dense(M, threshold)
        % CSC.from_dense() - Returns a structure containing a CSC
        % representation of the given dense matrix M
        % 
        % Construction of the sparse matrix is done by passing the transposed
        % of M to CSR.from_dense()
        %
        % INPUTS:
        %   - M: Matrix to be converted
        %   - threshold: Any number smaller than this value is considered 0.
        %                This in an optional argument that defaults to 0.
        %
        % OUTPUS:
        %   - M_sparse: Structure containing the sparse representation of M
        %       - val: Vector containing the values of the non-zero elements
        %       - row: Vector containing the row-indexes of each of the
        %              non-zero elements
        %       - col: Vector containing the column-index range of the
        %              non-zero elements
        %       - nnz: Number of non-zero elements
        %       - nrow: number of rows of M
        %       - ncol: number of columns of M
        if nargin == 1
            threshold = 0; % Any number smaller than threshold is considered a 0
        end

        % Compute dimensions of M
        [n, m] = size(M);
        
        %% Compute the CSR representation of the transposed
        [~, val, col, row, nnz] = CSR.from_dense(M', threshold); % NOTE: the col and row outputs are switched
        
        % Construct the output structure
        M_sparse.val = val;
        M_sparse.row = row;
        M_sparse.col = col;
        M_sparse.nnz = nnz;
        M_sparse.nrow = n;
        M_sparse.ncol = m;

    end

    end % Static methods

end
