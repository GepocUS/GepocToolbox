%% CSR - Class for representing and working with sparse matrices
% in the Compressed Sparse Row format
% 
% TODO:
%   - Add support for operation CSR*matrix
%   - Add support for operation vector*CSR
%   - Add support for operation matrix*CSR
%   - Add support for operation CSR*CSR
%   - Add support for operation CSR*CSR
%   - Add support for addition/subtraction operations (CSR+CSR, CSR+full, etc.)
%
% This class is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox

classdef CSR < GT_SparseMatrix

    methods

    %% CONSTRUCTOR

    function self = CSR(M, threshold)
        % CSR constructor
        % INPUTS:
        %   - M: Matrix to be converted to CSR format
        %   - threshold: Any number smaller than this value is considered 0.
        %                This in an optional argument that defaults to 0.
        if nargin == 0
            return;
        elseif nargin == 1
            threshold = 0; % Any number smaller than threshold is considered a 0
        end

        [~, self.val, self.row, self.col, self.nnz, self.nrow, self.ncol] = CSR.from_dense(M, threshold);

    end

    function M = to_dense(self)
        % CSR.to_dense - Return the dense representation of the CSR matrix
    
        % Initialize M
        M = zeros(self.nrow, self.ncol);

        % Fill dense matrix
        for i = 1:self.nrow
            for j = self.row(i):self.row(i+1)-1
                M(i, self.col(j)) = self.val(j);
            end
        end

    end

    %% OPERATORS

    function r = mtimes(obj, M)
        % CSR.mtimes - Operator overload for multiplication operations

        % Check that the multiplication is of self*M
        if ~isa(obj, 'CSR')
            error("CSR.mtimes: Premultiplications of CSR by a vector/matrix not currently supported.")
        end

        % Get size of M
        [n_M, m_M] = size(M);
        if m_M > 1
            error("CSR.mtimes: Currently only supports multiplication with column vector.");
        end
        if n_M ~= obj.ncol
            error("CSR.mtimes: Incorrect dimensions for matrix multiplication.");
        end

        % Initialize
        n = length(obj.row)-1;
        r = zeros(n, 1);
        
        % Perform multiplication
        for i = 1:n
            for j = obj.row(i):obj.row(i+1)-1
                    r(i) = r(i) + obj.val(j)*M(obj.col(j));
            end
        end

    end

    function r = ctranspose(obj)
        % CSR.mtranspose - Operator overload for transpose operator
        % TODO: Finish this method

        % Initialize
        r = CSR();
        r.nnz = obj.nnz;
        r.nrow = obj.ncol;
        r.ncol = obj.nrow;
        r.val = zeros(obj.nnz, 1);
        r.row = obj.ncol + 2;
        r.col = obj.nnz;

    end

    end % Methods

    methods (Static)

    function [M_sparse, val, row, col, nnz, n, m] = from_dense(M, threshold)
        % CSR.from_dense() - Returns a structure containing a CSR
        % representation of the given dense matrix M
        %
        % INPUTS:
        %   - M: Matrix to be converted
        %   - threshold: Any number smaller than this value is considered 0.
        %                This in an optional argument that defaults to 0.
        %
        % OUTPUS:
        %   - M_sparse: Structure containing the sparse representation of M
        %       - val: Vector containing the values of the non-zero elements
        %       - row: Vector containing the row-index range of the
        %              non-zero elements
        %       - col: Vector containing the column-indexes of each of the
        %              non-zero elements
        %       - nnz: Number of non-zero elements
        %       - nrow: number of rows of M
        %       - ncol: number of columns of M
        if nargin == 1
            threshold = 0; % Any number smaller than threshold is considered a 0
        end

        %% Initialize
        [n, m] = size(M); % Dimensions of M
        nnz = 0; % Number of non-zero elements
        val = [];
        col = [];
        row = zeros(1, n+1);
        
        %% Compute the CSR form
        for i = 1:n
            row(i) = nnz + 1;
            for j = 1:m
                if abs(M(i, j)) > threshold
                    val = [val M(i, j)];
                    col = [col j];
                    nnz = nnz + 1;
                end
            end
        end
        
        row(n + 1) = nnz + 1; % Add the last element to row
        
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
