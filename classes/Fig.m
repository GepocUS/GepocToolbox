%% Fig - A class for creating nice-looking figures
% 
% This class is intended to override Matlab's figure() so that
% figures are created with (opinionated) sensible default values
% that provide a nice-looking figure.
%
% The class is designed to be able to be easily switched for the
% figure() function and immediately get the benefits of the class,
% i.e., so that you can substitute
%   figure(1); plot(x, y, 'r:')
% for
%   Fig(1); plot(x, y, 'r:')
% and the resulting figure will be nicer looking.
%
% fig_var = Fig(num) - Fig constructor. Creates a new Fig object
%
% Constructor optional arguments
%   - clear_fig: boolean (true). Does clf() on the figure
%   - title: string. Title of the figure
%   - xlabel: string. xlabel of the figure
%   - ylabel: string. ylabel of the figure
%   - interpreter: string ('latex'). Interpreter used for writing text
%   - bg_color: string ('w'). Figure background color
%   - grig: bool (true). If true it sets the grid to on
%   - minorgrid: bool (false). If true it sets the minor grid to on
%   - hold: bool (true). It true it sets hold to on
%   - linewidth: scalar (1.5). Sets the default linewidth for the plots
%   - markersize: scalar (4). Sets the default markersize for the plots
%   - fontsize: scalar (20). Sets the default font size for text
%   - colorscheme: string ("lines"). Determines the color scheme of the plots
%   - max_num_colors: integer (8). Determines the number of different colors
%
% Fig properties (accessed using fig_var.property_name)
%   - All the contructor optional arguments are saved into 
%     properties with the same name except for: clear_fig
%   - num: Stores the number of the figure
%   - fh: Handler to the figure
%   - ax: Handler for the axis object of the figure
%   - ph: Cell containing the handlers to each of the plots
%
% Most properties can be reasigned: fig_var.property_name = value;
%
% Fig methods (executed using fig_var.method_name(args);)
%   - focus() focuses the figure
%   - clear() equivalent to clf(fig_var.num)
%   - trim() to delete (or select) plot margins
%   - y_scale() Switch between 'log' and 'linear' scale in y axis
%   - x_scale() Switch between 'log' and 'linear' scale in x axis
%   - plot() overrides Matlab's default plot() function
%   - save() saves the figure as a file
%
% For additional help on a method call: help Fig.method_name
%
% See also: figure, plot, clf

classdef Fig < handle
    
    properties
        title {ischar} % Title of the figure
        xlabel {ischar} = "" % Label of the X-axis
        ylabel {ischar} = "" % Label of the Y-axis
        bg_color % Background color
        grid {mustBeInteger, mustBeGreaterThanOrEqual(grid,0), mustBeLessThanOrEqual(grid,1)} % Sets grid to 'on' or 'off'
        hold {mustBeInteger, mustBeGreaterThanOrEqual(hold,0), mustBeLessThanOrEqual(hold,1)} % Sets hold to 'on' or 'off'
        minorgrid {mustBeInteger, mustBeGreaterThanOrEqual(minorgrid,0), mustBeLessThanOrEqual(minorgrid,1)} % Sets minor grid 'on' or 'off'
        linewidth {mustBeReal, mustBePositive} % Line width for new plot lines
        markersize {mustBeReal, mustBePositive} % Marker size for new plot lines
        fontsize {mustBeReal, mustBePositive} % Font size of the main text elements
        colorscheme {ischar} = "lines" % Color scheme of the figure
    end
    properties(SetAccess=protected, GetAccess=public)
       num % Stores the figure number
    end
    properties (Hidden = true)
        fh % Figure handler
        ax % Axis handler
        ph % List of plot handlers
        max_num_colors {mustBeInteger, mustBeGreaterThanOrEqual(max_num_colors,0)} = 8 % Number of plot lines before colors start repeating
        interpreter {ischar} = 'latex' % Interpreter used to print text
    end
    properties (Hidden = true, SetAccess=protected, GetAccess=protected)
        color_red = [1 0 0]; % Basic red color
        color_green = [0 1 0]; % Basic green color
        color_blue = [0 0 1]; % Basic blue color
        color_cyan = [0 1 1]; % Basic cyan color
        color_magenta = [1 0 1]; % Basic magenta color
        color_yellow = [1 1 0]; % Basic yellow color
        color_black = [0 0 0]; % Basic black color
        color_white = [1 1 1]; % Basic white color
        previous_ax_position; % Used to store the previous value of ax.Position. Used in previous_pos()
        default_date_Format = 'yy_MM_dd_HH_mm_ss';
    end
    properties (Hidden = true, SetAccess=protected)
        id % A unique id set when the object is created
        date % Instance of datetime created at the end of Fig's constructor
        gt_version % The version of GepocToolbox
        gt_id % The git hash version of GepocToolbox
    end
    
    methods

    %% CONSTRUCTOR
    
    function self = Fig(varargin)
        
        % Default values
        def_fig_num = [];
        def_clear_fig = true;
        def_title = "";
        def_xlabel = "";
        def_ylabel = "";
        def_interpreter = 'latex';
        def_bg_color = 'w';
        def_grid = true;
        def_hold = true;
        def_minorgrid = false;
        def_linewidth = 1.5;
        def_markersize = 4;
        def_fontsize = 20;
        def_max_num_colors = 8;
        def_colorscheme = "lines";
        
        % Parser
        par = inputParser;
        par.CaseSensitive = false;
        par.FunctionName = 'gpFig_constructor';
        
        % Optional
        addOptional(par, 'fig_num', def_fig_num, @(x) isnumeric(x) && (x>=1) && x==floor(x));
        % Name-value parameters
        addParameter(par, 'clear_fig', def_clear_fig, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'title', def_title, @(x) ischar(x));
        addParameter(par, 'xlabel', def_xlabel, @(x) ischar(x));
        addParameter(par, 'ylabel', def_ylabel, @(x) ischar(x));
        addParameter(par, 'interpreter', def_interpreter, @(x) ischar(x));
        addParameter(par, 'bg_color', def_bg_color, @(x) ischar(x) || isnumeric(x));
        addParameter(par, 'grid', def_grid, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'hold', def_hold, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'minorgrid', def_minorgrid, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'linewidth', def_linewidth, @(x) isnumeric(x) && (x>0));
        addParameter(par, 'markersize', def_markersize, @(x) isnumeric(x) && (x>0));
        addParameter(par, 'fontsize', def_fontsize, @(x) isnumeric(x) && (x>0) && x==floor(x));
        addParameter(par, 'max_num_colors', def_max_num_colors, @(x) mod(x,1)==0 && (x>0));
        addParameter(par, 'colorscheme', def_colorscheme, @(x) ischar(x));

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
        self.xlabel = par.Results.xlabel;
        self.ylabel = par.Results.ylabel;
        self.bg_color = par.Results.bg_color;
        self.grid = par.Results.grid;
        self.minorgrid = par.Results.minorgrid;
        self.linewidth = par.Results.linewidth;
        self.markersize = par.Results.markersize;
        self.fontsize = par.Results.fontsize;
        self.max_num_colors = par.Results.max_num_colors;
        self.colorscheme = par.Results.colorscheme;
        self.ph = cell(0);
        self.previous_ax_position = self.ax.Position;
        
        % Set identifiers and date
        self.id = string(java.util.UUID.randomUUID);
        try
            [self.gt_version, self.gt_id] = gepoc('version'); 
        catch ME
            if ~(strcmp(ME.identifier,'MATLAB:UndefinedFunction'))
                rethrow(ME)
            end
        end
        self.date = datetime('now');
        self.date.Format = self.default_date_Format;
        
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
    
    function set.xlabel(self, value)
        if ~isempty(value)
            self.xlabel = value;
            set(self.ax.XLabel, 'String', value);
        end
    end
    
    function set.ylabel(self, value)
        if ~isempty(value)
            self.ylabel = value;
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
    
    function set.linewidth(self, value)
        if ~isempty(value)
            self.linewidth = value;
            set(self.fh, 'DefaultLineLineWidth', value);
            self.update_plots_linewidth(value); % Update the line width of all plots
        end
    end
    
    function set.fontsize(self, value)
        if ~isempty(value)
            self.fontsize = value;
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
    
    function set.minorgrid(self, value)
        self.minorgrid = value;
        if value == true
            set(self.ax, 'XMinorGrid', 'on');
            set(self.ax, 'YMinorGrid', 'on');
        else
            set(self.ax, 'XMinorGrid', 'off');
            set(self.ax, 'YMinorGrid', 'off');
        end
    end

    function set.colorscheme(self, value)
        color_func = [value + "(self.max_num_colors)"];
        newColors = eval(color_func);
        self.colorscheme = value;
        self.ax.ColorOrder = newColors;
        self.update_plots_color(); % Update the color of all plots
    end

    function set_dateFormat(self, value)
        % This method sets the format of the internal date property
        % (instance of datetime) to the provided value
        self.date.Format = value;
    end
    
    %% PUBLIC METHODS
    
    function focus(self)
        % Fig.focus() - Focuses the figure
        %
        % Equivalent to calling figure(x) for some preexisting figure number x

        figure(self.fh);
    end
    
    function clear(self)
        % Fig.clear() - Clears the figure
        %
        % Calls clf() on the figure

        clf(self.num);
    end
    
    function trim(self, varargin)
        % Fig.trim() - Trims the empty space at the edges of the figure
        % 
        % Useful for making figures with no extra space for inserting them into articles
        % 
        % Fig.trim() trims the figure leaving no margin
        % 
        % Fig.trim(margin) trims the figure leaving the given margin
        % 
        % Fig.trim('margin_name', value) sets the provided margin to the given value
        % Possible margins are: 'margin' ('m') Same as Fig.trim(margin)
        %                       'west_margin' ('west', 'w') Left margin
        %                       'east_margin' ('east', 'e') Right margin 
        %                       'south_margin' ('south', 's') Bottom margin 
        %                       'north_margin' ('north', 'n') Top margin 
        % Specific margins take precedence over the value of 'margin'
        %
        % Fig.trim('undo') calls Fig.previous_pos(). Will undo the trim if called before
        % some other method which updates the plot position (see Fig.previous_pos())
        %
        % See also: Fig.previous_pos

        if nargin == 2 && strcmp(varargin{1}, 'undo')

            self.previous_pos();
            return;

        elseif nargin == 2 && isnumeric(varargin{1})

            margin = varargin{1};
            west_margin = margin;
            east_margin = margin;
            south_margin = margin;
            north_margin = margin;

        else

            % Default values
            def_margin = 0.0;
            def_west_margin = NaN;
            def_east_margin = NaN;
            def_south_margin = NaN;
            def_north_margin = NaN;

            % Parser
            par = inputParser;
            par.CaseSensitive = false;
            par.FunctionName = 'Fig.trim()';
            % Name-value parameters
            addParameter(par, 'margin', def_margin, @(x) isnumeric(x));
            addParameter(par, 'west_margin', def_west_margin, @(x) isnumeric(x));
            addParameter(par, 'east_margin', def_east_margin, @(x) isnumeric(x));
            addParameter(par, 'south_margin', def_south_margin, @(x) isnumeric(x));
            addParameter(par, 'north_margin', def_north_margin, @(x) isnumeric(x));
            % Parse
            parse(par, varargin{:})
            % Rename and set
            margin = par.Results.margin;
            west_margin = margin;
            east_margin = margin;
            south_margin = margin;
            north_margin = margin;
            if ~isnan(par.Results.west_margin);  west_margin = par.Results.west_margin; end
            if ~isnan(par.Results.east_margin);  east_margin = par.Results.east_margin; end
            if ~isnan(par.Results.south_margin); south_margin = par.Results.south_margin; end;
            if ~isnan(par.Results.north_margin); north_margin = par.Results.north_margin; end;

        end

        outerpos = self.ax.OuterPosition;
        ti = self.ax.TightInset; 
        left = outerpos(1) + ti(1) + west_margin;
        bottom = outerpos(2) + ti(2) + south_margin;
        ax_width = outerpos(3) - ti(1) - ti(3) - east_margin - west_margin;
        ax_height = outerpos(4) - ti(2) - ti(4) - north_margin - south_margin;
        self.ax.Position = [left bottom ax_width ax_height];

    end

    function save(self, varargin)
        % Fig.save() - Saves the figure
        %
        % This method saves the figure with the given name and
        % file-type provided in the provided directory.
        %
        % Name-value input parameters:
        %   - name: Name of the file. If none is provided then the figure's title is used.
        %           If no title is available then the date of creation of the Fig is used.
        %   - directory: Directory where the file is saved.
        %                Relative paths can be used by "./rel_path" or simply "rel_path"
        %                Warning: "~" cannot be used in Linux. Use "/home/username/"
        %   - extension: Format the file is saved as. By default it is saved as a '.eps'
        %                Supported file-types are the ones supported by saveas()
        %   - date: If true, the date of creation of the figure is appended to the name
        %   - unique: If true, the Fig unique identifier is appended to the name
        %   - color: defaults to true. If false then the figure is saved in black and white.
        %            Only available for the 'eps' extension.
        %
        % See also: saveas

        % Default values
        def_name = [];
        def_directory = './';
        def_extension = 'eps';
        def_date = false;
        def_unique = false;
        def_color = true;

        % Parser
        par = inputParser;
        par.CaseSensitive = false;
        par.FunctionName = 'Fig.save()';
        % Name-value parameters
        addParameter(par, 'name', def_name, @(x) isstring(x) || ischar(x));
        addParameter(par, 'directory', def_directory, @(x) isstring(x) || ischar(x));
        addParameter(par, 'extension', def_extension, @(x) ischar(x));
        addParameter(par, 'date', def_date, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'unique', def_unique, @(x) islogical(x) || x==1 || x==0);
        addParameter(par, 'color', def_color, @(x) islogical(x) || x==1 || x==0);
        % Parse
        parse(par, varargin{:})

        % Set name
        name = string(par.Results.name);
        if isempty(name)
            if ~isempty(self.title)
                name = self.title;
                if par.Results.date
                    name = name + "_" + string(self.date);
                end
            else
                name = string(self.date);
            end
        else
            if par.Results.date
                name = name + "_" + string(self.date);
            end
        end

        % Add unique identifier
        if par.Results.unique
            name = name + "_" + extractBetween(self.id, 1, "-");
        end

        % Create full path
        directory = string(par.Results.directory);
        dir_char = char(directory);
        % directory
        if ~strcmp(dir_char(end), "/")
            directory = directory + "/";
        end

        name = fullfile(directory, name);

        % Save figure
        switch par.Results.extension
            case 'eps'
                if par.Results.color
                    saveas(self.fh, name, 'epsc');
                else
                    saveas(self.fh, name, 'eps');
                end
            otherwise
                saveas(self.fh, name, par.Results.extension);
        end

    end

    function y_scale(self, value)
        % Fig.y_scale() - Select y axis scale mode
        % 
        % y_scale() - Switch Y axis scale between 'linear' and 'log'
        % y_scale('log') - Set Y axis scale to 'log'
        % y_scale('linear') - Set Y axis scale to 'linear'
        %
        % See also: Fig.x_scale

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
        % Fig.x_scale() - Select x axis scale mode
        %
        % x_scale() - Switch X axis scale between 'linear' and 'log'
        % x_scale('log') - Set X axis scale to 'log'
        % x_scale('linear') - Set X axis scale to 'linear'
        %
        % See also: Fig.y_scale

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
        % Fig.plot() - Overload of Matlab's standard plot() function
        %
        % Supports the main name-value parameters of Matlab's plot(), but
        % providing (opinionated) sensible values to them.
        % 
        % The currently supported name-value parameters are:
        %   - linewidth
        %   - markersize
        %   - linestyle
        %   - marker
        %   - color
        % 
        % Also supports the standard Fig.plot(x, y, 'r:') way of
        % choosing line color and style.
        %
        % See also: plot

        % Default values
        def_mods = '';
        def_linewidth = self.linewidth;
        def_markersize = self.markersize;
        def_linestyle = '-';
        def_marker = 'none';
        def_color = self.ax.ColorOrder(mod(length(self.ph), self.max_num_colors) + 1, :);

        % Parser
        par = inputParser;
        par.CaseSensitive = false;
        par.FunctionName = 'Fig.plot()';
        % Required
        addRequired(par, 'x');
        addRequired(par, 'y');
        % Optional
        addOptional(par, 'mods', def_mods, @(x) ischar(x));
        % Name-value parameters
        addParameter(par, 'linestyle', def_linestyle, @(x) ischar(x));
        addParameter(par, 'marker', def_marker, @(x) ischar(x));
        addParameter(par, 'linewidth', def_linewidth, @(x) isnumeric(x) && (x>0));
        addParameter(par, 'markersize', def_markersize, @(x) isnumeric(x) && (x>0));
        addParameter(par, 'color', def_color);
        % Parse
        if mod(length(varargin), 2)==0
            parse(par, varargin{1}, varargin{2}, def_mods, varargin{3:end});
        else
            parse(par, varargin{:});
        end
        % Rename
        x = par.Results.x;
        y = par.Results.y;
        mods_ = par.Results.mods;
        linestyle_ = par.Results.linestyle;
        marker_ = par.Results.marker;
        linewidth_ = par.Results.linewidth;
        markersize_ = par.Results.markersize;
        color_ = par.Results.color;

        % Use marker from mods is available
        idx_mods_marker = regexp(mods_ ,'[.ox+*sdv^<>ph]');
        if ~isempty(idx_mods_marker)
            marker_ = mods_(idx_mods_marker);
            linestyle_ = 'none';
        end

        % Use color from mods if available
        idx_mods_color = regexp(mods_ ,'[rgbcmykw]');
        if ~isempty(idx_mods_color)
            color_ = self.get_basic_color(mods_(idx_mods_color));
        end

        % Use linestyle from mods if available
        mods_linestyle = erase(mods_, mods_(idx_mods_marker));
        mods_linestyle = erase(mods_linestyle, mods_(idx_mods_color));
        if ~isempty(mods_linestyle)
            linestyle_ = mods_linestyle;
        end

        % Plot
        self.focus();
        self.ph{end+1} = plot(self.ax, x, y, mods_,...
                              'Color', color_,...
                              'linewidth', linewidth_, 'LineStyle', linestyle_,...
                              'markersize', markersize_, 'Marker', marker_...
                              );
         
        % Post plot
        self.previous_ax_position = self.ax.Position;

    end
    
    end % End public methods

    %% PROTECTED METHODS

    methods (Access = protected)

    function previous_pos(self)
        % Fig.previous_pos() - Recover previous axis position
        % 
        % Returns the figure axes position to its previous stored value
        % Functions that update the previous value are:
        %   Fig.plot(), Fig.trim()
        %
        % See also: Fig.plot, Fig.trim

        pos_aux = self.ax.Position;
        self.ax.Position = self.previous_ax_position;
        self.previous_ax_position = pos_aux;
    end

    function update_plots_linewidth(self, value)
        % Sets the linewidth of all the current plots to the given value
        for i = 1:length(self.ph)
            self.ph{i}.LineWidth = value;
        end
    end

    function update_plots_color(self)
        % Sets the color order of all the current plots to the colors in ax.ColorOrder
        for i = 1:length(self.ph)
            self.ph{i}.Color = self.ax.ColorOrder(mod(i-1, self.max_num_colors)+1, :);
        end
    end

    function color = get_basic_color(self, name)
        % Returns the basic color states in the string 'name'
        color = [];

        switch name
            case 'red';     color = self.color_red;
            case 'r';       color = self.color_red;
            case 'green';   color = self.color_green;
            case 'g';       color = self.color_green;
            case 'blue';    color = self.color_blue;
            case 'b';       color = self.color_blue;
            case 'cyan';    color = self.color_cyan;
            case 'c';       color = self.color_cyan;
            case 'magenta'; color = self.color_magenta;
            case 'm';       color = self.color_magenta;
            case 'yellow';  color = self.color_yellow;
            case 'y';       color = self.color_yellow;
            case 'black';   color = self.color_black;
            case 'k';       color = self.color_black;
            case 'white';   color = self.color_white;
            case 'w';       color = self.color_white;
        end

    end

    end % End of protected methods
    
end
