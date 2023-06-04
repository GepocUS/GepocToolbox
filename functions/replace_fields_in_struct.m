%% replace_fields_in_struct
%
% original = replace_fields_in_struct(original, struct_in)
% 
% Replaces the fields in original that by the same fields in struct_in.
% Fields that are in original but not in struct_in are unchanged.
% Fields that are in struct_in but not in original are ignored.
% 
% It returns the 'original' structure with the shared fields replaced.
% 
% This function is part of GepocToolbox: https://github.com/GepocUS/GepocToolbox
%

function original = replace_fields_in_struct(original, struct_in)

    fn = fieldnames(original); % Get the field names
    
    % Override fields
    for i = 1:numel(fn)
        if isfield(struct_in, fn{i})
            original.(fn{i}) = struct_in.(fn{i});
        end
    end

end