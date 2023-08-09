%% gt_exec_tests - Execute the GepocToolbox unit tests
%
% This function executes a series of "unit" tests to determine that the
% functions and classes of the toolbox are working correctly.

function out = gt_exec_tests(opt)

    % Check that the toolbox is installed by searching one of its functions
    if ~exist('compute_PIS_ellip.m')
        warning("GepocToolbox not correctly installed. Run: gepoc('install')");
        return
    end

    % Default opt
    def_opt.verbose = 1;

    if nargin == 0
        opt = def_opt;
    end
    opt = add_default_options_to_struct(opt, def_opt);

    % Get list of unit tests
    i = 1;

    test{i}.fcn = "gt_test_sparse";
    test{i}.desc = "Test functions related to sparse matrix representation and operations";
    i = i+1;

    % Run unit tests
    tt = cputime;
    n_tests = length(test);

    for i = 1:n_tests
        try
            t=cputime;           
            if opt.verbose
                fprintf("Testing function: " + test{i}.fcn + ". ");
            end
            [pass(i), info(i)] = eval(test{i}.fcn + "(opt)");
            ttime(i) = cputime-t;
            if opt.verbose
                if pass(i)
                    fprintf("Results: pass\n");
                else
                    fprintf("Results: FAIL\n");
                end
            end
        catch e
            pass(i) = false;   
            info{i} = struct('pass', false, 'msg', e.message);
            ttime(i) = cputime-tt;
            fprintf("Results: FAIL\n");
        end
        % Print messages
        if opt.verbose > 0
            for j = 1:length(info(i).pass)
                fprintf("\tSubtest: " + info(i).desc(j) + ". ");
                if info(i).pass(j) == true
                    fprintf("Results: pass. \n");
                else
                    fprintf("Results: FAIL. \n");
                end
                fprintf("\t\t Msg: " + info(i).msg{j} + "\n");
            end
        end
    end
    totaltime = cputime-tt;

    % Analyze results
    num_pass = sum(pass);
    num_err = n_tests - num_pass;

    if num_err == 0
        fprintf("\nAll " + num2str(n_tests) + " tests PASSED!\n");
    else
        fprintf("\n" + num2str(num_err) + " tests FAILED out of " + num2str(n_tests) + " tests.\n");
        fprintf("\nThe following tests failed:\n");
        for i = 1:n_tests
            if pass(i) == false
                fprintf(test{i}.fcn + "\n");
                for j = 1:length(info(i).pass)
                    if info(i).pass(j) == false
                        fprintf("\t" + info(i).msg{j} + "\n");
                    end
                end
            end
        end
    end

    % Print error messages
    if num_err~=0
    end

    % Return information
    out.pass = pass;
    out.info = info;

end
