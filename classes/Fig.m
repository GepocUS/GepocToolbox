%% A class for creating nice-looking figures
% 
% Constructor optional arguments
%   - clear_fig: boolean (true). Does clf() on the figure
%   - title: string. Title of the figure
%   - x_label: string. x_label of the figure
%   - y_label: string. y_label of the figure
%   - interpreter: string ('latex'). Interpreter used for writing text
%   - bg_color: string ('w'). Figure background color
%   - grig: bool (true). If true it sets the grid to on
%   - minor_grid: bool (false). If true it sets the minor grid to on
%   - hold: bool (true). It true it sets hold to on
%   - line_width: scalar (1.5). Sets the default line_width for plots
%   - font_size: scalar (20). Sets the default font size for text
%   - color_scheme: string ("lines"). Determines the color scheme of the plots
%   - max_num_colors: integer (8). Determines the number of different colors
%
% Properties
%   - All the contructor optional arguments are saved into 
%     properties with the same name except for: clear_fig
%   - num: Stores the number of the figure
%   - fh: Handler to the figure
%   - ax: Handler for the axis object of the figure
%   - ph: Cell containing the handlers to each of the plots

classdef Fig < handle
    
    properties
        title {ischar} % Title of the figure
        x_label {ischar} = "" % Label of the X-axis
        y_label {ischar} = "" % Label of the Y-axis
        bg_color % Background color
        grid {mustBeInteger, mustBeGreaterThanOrEqual(grid,0), mustBeLessThanOrEqual(grid,1)} % Sets grid to 'on' or 'off'
        hold {mustBeInteger, mustBeGreaterThanOrEqual(hold,0), mustBeLessThanOrEqual(hold,1)} % Sets hold to 'on' or 'off'
        minor_grid {mustBeInteger, mustBeGreaterThanOrEqual(minor_grid,0), mustBeLessThanOrEqual(minor_grid,1)} % Sets minor grid 'on' or 'off'
        line_width {mustBeReal, mustBePositive} % Line width for new plot lines
        font_size {mustBeReal, mustBePositive} % Font size of the main text elements
        color_scheme {ischar} = "lines"
    end
    properties(SetAccess=protected, GetAccess=public)
       num % Stores the figure number
    end
    properties (Hidden = true)
        fh % Figure handler
        ax % Axis handler
        ph % List of plot handlers
        max_num_colors {mustBeInteger, mustBeGreaterThanOrEqual(max_num_colors,0)} = 8
        interpreter {ischar} = 'latex' % Interpreter used to print text
    end
    
    methods

    %% CONSTRUCTOR
    
    function self = Fig(varargin)
        
        % Default values
        def_fig_num = [];
        def_clear_fig = true;
        def_title = "";
        def_x_label = "";
        def_y_label = "";
        def_interpreter = 'latex';
        def_bg_color = 'w';
        def_grid = true;
        def_hold = true;
        def_minor_grid = false;
        def_line_width = 1.5;
        def_font_size = 20;
        def_max_num_colors = 8;
        def_color_scheme = "lines";
        
        % Parser
        par = inputParser;
        par.CaseSensitive = false;
        par.FunctionName = 'gpFig_constructor';
        
        % Optional
        addOptional(par, 'fig_num', def_fig_num, @(x) isnumeric(x) && (x>=1) && x==floor(x));
        % Name-value parameters
        addParameter(par, 'clear_fig', def_clear_fig, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'title', def_title, @(x) ischar(x));
        addParameter(par, 'x_label', def_x_label, @(x) ischar(x));
        addParameter(par, 'y_label', def_y_label, @(x) ischar(x));
        addParameter(par, 'interpreter', def_interpreter, @(x) ischar(x));
        addParameter(par, 'bg_color', def_bg_color, @(x) ischar(x) || isnumeric(x));
        addParameter(par, 'grid', def_grid, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'hold', def_hold, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'minor_grid', def_minor_grid, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'line_width', def_line_width, @(x) isnumeric(x) && (x>0));
        addParameter(par, 'font_size', def_font_size, @(x) isnumeric(x) && (x>0) && x==floor(x));
        addParameter(par, 'max_num_colors', def_max_num_colors, @(x) mod(x,1)==0 && (x>0));
        addParameter(par, 'color_scheme', def_color_scheme, @(x) ischar(x));

        % Parse
        parse(par, varargin{:})
        % Rename
        fig_num = par.Results.fig_num;
        clear_fig = par.Results.clear_fig;
        
        % Create figure
        if isempty(fig_num)
            self.fh = figure();
        else
            self.fh = figure(fig_num);
        end
        if clear_fig; clf(self.fh.Number); end % Clear figure
        
        self.ax = gca; % Get axis handler
        
        % Set properties
        self.hold = par.Results.hold;
        self.interpreter = par.Results.interpreter;
        self.title = par.Results.title;
        self.x_label = par.Results.x_label;
        self.y_label = par.Results.y_label;
        self.bg_color = par.Results.bg_color;
        self.grid = par.Results.grid;
        self.minor_grid = par.Results.minor_grid;
        self.line_width = par.Results.line_width;
        self.font_size = par.Results.font_size;
        self.max_num_colors = par.Results.max_num_colors;
        self.color_scheme = par.Results.color_scheme;
        self.ph = cell(0);
        
    end
    
    %% GETTERS and SETTERS
    
    function value = get.num(self)
        value = self.fh.Number;
    end
    
    function set.title(self, value)
        if ~isempty(value)
            self.title = value;
            set(self.ax.Title, 'String', value);
        end
    end
    
    function set.x_label(self, value)
        if ~isempty(value)
            self.x_label = value;
            set(self.ax.XLabel, 'String', value);
        end
    end
    
    function set.y_label(self, value)
        if ~isempty(value)
            self.y_label = value;
            set(self.ax.YLabel, 'String', value);
        end
    end
    
    function set.interpreter(self, value)
        self.interpreter = value;
        set(self.ax.Title,'Interpreter', value);
        set(self.ax, 'TickLabelInterpreter', value);
        set(self.ax.XLabel, 'Interpreter', value);
        set(self.ax.YLabel, 'Interpreter', value);
    end
    
    function set.line_width(self, value)
        if ~isempty(value)
            self.line_width = value;
            set(self.fh, 'DefaultLineLineWidth', value);
            self.update_plots_line_width(value); % Update the line width of all plots
        end
    end
    
    function set.font_size(self, value)
        if ~isempty(value)
            self.font_size = value;
            set(self.ax, 'FontSize', value);
        end
    end
    
    function set.bg_color(self, value)
        if ~isempty(value)
            self.bg_color = value;
            set(self.fh, 'Color', value);
        end
    end
    
    function set.grid(self, value)
        self.grid = value;
        if value == true
            set(self.ax, 'XGrid', 'on');
            set(self.ax, 'YGrid', 'on');
        else
            set(self.ax, 'XGrid', 'off');
            set(self.ax, 'YGrid', 'off');
        end
    end
    
    function set.hold(self, value)
        self.hold = value;
        if value == true
            hold(self.ax, 'on');
        else
            hold(self.ax, 'off');
        end
    end
    
    function set.minor_grid(self, value)
        self.minor_grid = value;
        if value == true
            set(self.ax, 'XMinorGrid', 'on');
            set(self.ax, 'YMinorGrid', 'on');
        else
            set(self.ax, 'XMinorGrid', 'off');
            set(self.ax, 'YMinorGrid', 'off');
        end
    end

    function set.color_scheme(self, value)
        color_func = [value + "(self.max_num_colors)"];
        newColors = eval(color_func);
        self.color_scheme = value;
        self.ax.ColorOrder = newColors;
        self.update_plots_color(); % Update the color of all plots
    end

    
    %% PUBLIC METHODS
    
    function focus(self)
        % Fig.focus() - Focuses the figure
        % Equivalent to calling figure(x) for some preexisting figure number x
        figure(self.fh);
    end
    
    function clear(self)
        % Fig.clear() - Clears the figure
        % Calles clf() on the figure
        clf(self.num);
    end
    
    function trim(self)
        % Fig.trim() - Trims the empty space at the edges of the figure
        % Useful for making figures with no extra space for inserting them into articles
        outerpos = self.ax.OuterPosition;
        ti = self.ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        self.ax.Position = [left bottom ax_width ax_height];
    end
    
    function y_scale(self, value)
        % y_scale() - Switch Y axis scale between 'linear' and 'log'
        % y_scale('log') - Set Y axis scale to 'log'
        % y_scale('log') - Set Y axis scale to 'linear'
        if nargin == 1
            if strcmp(self.ax.YScale, 'linear')
                value = 'log';
            else
                value = 'linear';
            end
        end
        set(self.ax, 'YScale', value);
    end
    
    function x_scale(self, value)
        % x_scale() - Switch X axis scale between 'linear' and 'log'
        % x_scale('log') - Set X axis scale to 'log'
        % x_scale('log') - Set X axis scale to 'linear'
        if nargin == 1
            if strcmp(self.ax.XScale, 'linear')
                value = 'log';
            else
                value = 'linear';
            end
        end
        set(self.ax, 'XScale', value);
    end
    
    
    %% PLOT METHODS
    
    function plot(self, varargin)
        % Overload of the standard plot() function
        % TODO: work on default values, etc.
        % TODO: seamless integration of markers

        % Default values
        def_linewidth = self.line_width;
        def_markersize = 4;
        def_color = [];
        


        self.focus();
        self.ph{end+1} = plot(self.ax, varargin{:}, 'linewidth', self.line_width, 'markersize', 4);
    end
    
    end % End public methods

    %% PROTECTED METHODS

    methods (Access = protected)

    function update_plots_line_width(self, value)
        for i = 1:length(self.ph)
            self.ph{i}.LineWidth = value;
        end
    end

    function update_plots_color(self)
        for i = 1:length(self.ph)
            self.ph{i}.Color = self.ax.ColorOrder(mod(i-1, self.max_num_colors)+1, :);
        end
    end

    end
    
end
