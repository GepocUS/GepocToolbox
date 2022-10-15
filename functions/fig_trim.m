%% fig_trim - Trim empty space surrounding figure
% 
% This script runs a series of commands that trim the 
% empty space surrounding the last figure focused on
%
% This script is part of the GepocToolbox: https://github.com/GepocUS/GepocToolbox
% 

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
