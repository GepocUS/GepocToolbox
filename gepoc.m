%% gepoc - root function of the GepocToolbox toolbox
%
% This toolbox contains an assortment of classes and function useful
% for those working in the fields of optimization and predictive control.
% 
% For information on the toolbox please visit https://github.com/GepocUS/GepocToolbox
%
% This function can be used for various purposes:
%
%   gepoc('version')        % Returns the version number of the toolbox. If git is
%                             installed it also returns the git hash in the second
%                             output, i.e., [version, hash] = gepoc('version').
%                             The first output is the latest version tag.
%
%   gepoc('install')        % Installs the toolbox.
%   gepoc('uninstall')      % Uninstalls the toolbox.
%
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function varargout = gepoc(varargin)

if nargin == 0
    help gepoc
    return
end

switch varargin{1}
    
    
    % [version, git_hash] = gepoc('version')
    %
    % Returns information of the version of the toolbox.
    % 
    % OUTPUTS:
    %   - version: Version number of the toolbox -> v.MAJOR.MINOR.PATCH
    %   - git_hash: Current git hash of the toolbox. Only works if git
    %               is installed and the toolbox is being tracked.
    % 
    case 'version'
        
        varargout{1} = 'v0.0.1';
        
        % If git is installed it will return the hash of the current commit
        try
            [system_status, git_hash] = system('git rev-parse HEAD');
        catch
            varargout{2} = 'Could not obtain git hash';
        end
        
        if system_status == 0
            varargout{2} = convertCharsToStrings(git_hash(1:end-1));
        else
            varargout{2} = 'Could not obtain git hash';
        end
        
    % gepoc('install')
    %
    % Installs the toolbox. It adds the required directories to the
    % Matlab path.
    %
    case 'install'
        root_path = gt_get_root_directory();
        addpath(root_path);
        addpath([root_path '/classes/']);
        addpath([root_path '/functions/']);
        addpath([root_path '/solvers/']);
        addpath([root_path '/benchmarks/']);
        warning('off','MATLAB:rmpath:DirNotFound')
        addpath([root_path '/personal/']);
        warning('on','MATLAB:rmpath:DirNotFound')
        savepath
        disp('GepoxToolbox installed');
    
    % gepoc('uninstall')
    %
    % Uninstalls the toolbox. It removes the required directories from
    % the Matlab path.
    %
    case 'uninstall'
        root_path = gt_get_root_directory();
        warning('off','MATLAB:rmpath:DirNotFound')
        rmpath(root_path);
        rmpath([root_path '/classes/']);
        rmpath([root_path '/functions/']);
        rmpath([root_path '/solvers/']);
        rmpath([root_path '/benchmarks/']);
        warning('on','MATLAB:rmpath:DirNotFound')
        savepath
        disp('GepoxToolbox uninstalled');
    
    % Command not recognized
    otherwise
        varargout{1} = 'Command not recognized'; 
        
end

end
