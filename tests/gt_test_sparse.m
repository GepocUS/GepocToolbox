function [pass, info] = gt_test_sparse(opt)

    % Initialize
    pass = true;
    tt = 1; % Test counter

    % Construct random sparse matrix (in dense form)
    n = 6;
    m = 8;
    M = rand(n, m);
    M = M.*(M < 0.3); % Make the matrix sparse
    M(1, 2) = 0.7; % Make sure we have at least one non-zero element

    % Substest 1: Construct CSR and CSC instances
    
    info.desc{tt} = "Construct CSR matrix";
    try
        Mcsr = CSR(M);
        info.pass(tt) = true;
        info.msg{tt} = "CSR matrix constructed correctly";
    catch e
        info.pass(tt) = false;
        info.msg{tt} = "Error constructing CSR matrix: " + e.message;
    end
    tt = tt+1;

    info.desc{tt} = "Construct CSC matrix";
    try
        Mcsc = CSC(M);
        info.pass(tt) = true;
        info.msg{tt} = "CSC matrix constructed correctly";
    catch e
        info.pass(tt) = false;
        info.msg{tt} = "Error constructing CSC matrix: " + e.message;
    end
    tt = tt+1;
    
    % Subtest 2: Transform back to dense
    
    info.desc{tt} = "CSR to dense format";
    try
        Md = Mcsr.to_dense();
        if norm(Md - M, Inf) < eps
            info.pass(tt) = true;
            info.msg{tt} = "CSR conversion to dense format successful";
        else
            info.pass(tt) = false;
            info.msg{tt} = "Wrong CSR conversion to dense format";
        end
    catch e
        info.pass(tt) = false;
        info.msg{tt} = "Error converting CSR matrix to dense format: " + e.message;
    end
    tt = tt+1;
    
    % Subtest 3: Multiply CSR by a vector
    v = randn(m, 1);

    info.desc{tt} = "CSR*vector operation";
    try
        r = Mcsr*v;
        if norm(M*v - r, Inf) < eps
            info.pass(tt) = true;
            info.msg{tt} = "CSR*vector successful";
        else
            info.pass(tt) = false;
            info.msg{tt} = "Wrong result of CSR*vector";
        end
    catch e
        info.pass(tt) = false;
        info.msg{tt} = "Error in CSR*vector: " + e.message;
    end
    tt = tt+1;

    % Return results
    if any(info.pass == false)
        pass = false;
    end

end
