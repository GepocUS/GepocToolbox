%% gt_get_root_directory
%
% This function returns the root directory of the GepocToolbox.
% 
% This function is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

function gepoc_path = gt_get_root_directory()

    % Get the path to the directory where this function is saved
    full_path = mfilename('fullpath');
    gepoc_path = fileparts(full_path);
    
end
