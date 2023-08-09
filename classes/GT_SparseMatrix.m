%% GT_SparseMatrix - Abstract class for the CSC and CSR sparse matrix classes

classdef (Abstract) GT_SparseMatrix

    properties
        val % Non-zero values
        row % Row index information
        col % Column index information
        nnz % Number of non-zero elements
        nrow % Number of rows
        ncol % Number of columns
    end

    methods (Static, Abstract)
        from_dense(); % Returns the sparse matrix from the given dense matrix
    end

end
