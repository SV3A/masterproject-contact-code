classdef Plottools < handle
  % 'Plottools' is a class which defines a common plotting module.

  properties (Access = private)
    dbugplots = {};  % Cell array containing debug-plot objects

    orbit_plt % Orbit plot figure handle
    state_plt % State plot figure handle
    lat_plt   % Lateral plot figure handle
    int_plt   % Indentation plot figure handle
    intf_plt  % Indentation-force plot figure handle

    fontsize_lab = 28;
    fontsize_ax  = 26;
  end


  methods (Access = public)

    function obj = Plottools()
      % Constructor function
    end


    function debugplot(obj, varargin)
      % A public function handle for creating a debug plot.  The function takes
      % either 1 input (filename) or 8 (see 'debugplt1_internal').

      % If only one arg is given it should be a file to read inputs from
      if length(varargin) == 1
        obj.debugplt1_internal_f(varargin{1})
      elseif length(varargin) == 8
        obj.debugplt1_internal(varargin{:})
      else
        error('Too many arguments given')
      end
    end


    function orbit(obj, rot_x, rot_y, clearance, sup_x, sup_y)
      % Plots the rotor orbit inside the stator boundary (clearance circle).

      % Create figure or overwrite it
      obj.orbit_plt = obj.set_fig(obj.orbit_plt);

      % Set figure properties
      set(obj.orbit_plt, 'name', 'Orbit Plot', 'color', 'w', ...
          'units', 'normalized');

      % If MATLAB default win style is "docked" don't resize the figure
      if ~strcmp(get(0, 'DefaultFigureWindowStyle'), 'docked')
        set(obj.orbit_plt, 'outerposition', [0.5 0 0.5 1]);
      end

      % Rotor orbit (converted to [mm])
      plot(rot_x * 1e3, rot_y * 1e3, 'b', 'LineWidth', 1.2); grid on; hold on

      if nargin > 4
        plot(sup_x * 1e3, sup_y * 1e3, 'k', 'LineWidth', 1.2);
      end

      % THIS IS HACK
      %clearance = 0.00121;
      %clearance = 0.001106;
      clearance = 0.00114;

      % Clearance circle
      plot(clearance * cos(linspace(0, 2*pi)) * 1e3, ...
           clearance * sin(linspace(0, 2*pi)) * 1e3, 'r', 'LineWidth', 1.5)
      xlabel('Orbit [mm]', 'interpreter', 'latex', 'FontSize', 18)
      axis equal

      set(findall(gcf, '-property', 'FontSize'), ...
          'FontSize', obj.fontsize_ax, 'FontName', 'Times')
    end


    function lateral(obj, t, rot_x, rot_y)
      % Plot lateral displacement of the rotor

      % Create figure or overwrite it
      obj.lat_plt = obj.set_fig(obj.lat_plt);

      % Set figure properties
      set(obj.lat_plt, 'name', 'Lateral displacements', 'color', 'w', ...
          'units', 'normalized');

      if ~strcmp(get(0, 'DefaultFigureWindowStyle'), 'docked')
        set(obj.state_plt, 'outerposition', [0.5 0 0.5 1]);
      end

      subplot(2, 1, 1)
      plot(t, rot_x * 1e3, 'b', 'LineWidth', 1.2); grid on
      ylabel('$x$-displacement [mm]', 'interpreter', 'latex', 'FontSize', 18)

      subplot(2, 1, 2)
      plot(t, rot_y * 1e3, 'r', 'LineWidth', 1.2); grid on
      xlabel('Time [s]', 'interpreter', 'latex', ...
             'FontSize', obj.fontsize_lab)
      ylabel('$y$-displacement [mm]', 'interpreter', 'latex', ...
             'FontSize', obj.fontsize_lab)

      set(findall(gcf, '-property', 'FontSize'), ...
          'FontSize', obj.fontsize_ax, 'FontName', 'Times')
    end


    function indent(obj, time, delta, event_times)
      % Plot indentation displacement of the rotor

      if isempty(event_times)
        warning('WARNING No impacts to plot')
        return
      end

      % Create figure or overwrite it
      obj.int_plt = obj.set_fig(obj.int_plt);

      % Set figure properties
      set(obj.int_plt, 'name', 'Lateral displacements', 'color', 'w', ...
          'units', 'normalized');

      if ~strcmp(get(0, 'DefaultFigureWindowStyle'), 'docked')
        set(obj.state_plt, 'outerposition', [0.5 0 0.5 1]);
      end

      idx = zeros(length(event_times), 1);
      for i = 1:length(event_times)
        idx(i) = find(event_times(i) == time);
      end

      hold on

      for i = 2:2:length(event_times)
        plot(time(idx(i-1):idx(i)), delta(idx(i-1):idx(i))*1e6); grid on
      end
      xlabel('Time [s]', 'interpreter', 'latex', ...
             'FontSize', obj.fontsize_lab)
      ylabel('Indentation [$\mu$m]', 'interpreter', 'latex', ...
             'FontSize', obj.fontsize_lab)

      set(findall(gcf, '-property', 'FontSize'), ...
          'FontSize', obj.fontsize_ax, 'FontName', 'Times')
    end


    function indentf(obj, time, delta, fn, event_times)
      % Plot lateral displacement of the rotor

      if isempty(event_times)
        warning('WARNING No impacts to plot')
        return
      end

      % Create figure or overwrite it
      obj.intf_plt = obj.set_fig(obj.intf_plt);

      % Set figure properties
      set(obj.intf_plt, 'name', 'Lateral displacements', 'color', 'w', ...
          'units', 'normalized');

      if ~strcmp(get(0, 'DefaultFigureWindowStyle'), 'docked')
        set(obj.state_plt, 'outerposition', [0.5 0 0.5 1]);
      end

      idx = zeros(length(event_times), 1);
      for i = 1:length(event_times)
        idx(i) = find(event_times(i) == time);
      end

      hold on

      for i = 2:2:length(event_times)
        plot(delta(idx(i-1):idx(i))*1e6, fn(idx(i-1):idx(i))); grid on
      end
      xlabel('Indentation [$\mu$m]', 'interpreter', 'latex', ...
             'FontSize', obj.fontsize_lab)
      ylabel('Radial force [N]', 'interpreter', 'latex', ...
             'FontSize', obj.fontsize_lab)

      set(findall(gcf, '-property', 'FontSize'), ...
          'FontSize', obj.fontsize_ax, 'FontName', 'Times')
    end


    function states(obj, t, state_vector, y_labels)
      % Plots the state vector.
      % INPUT:
      %   t           : Time vector
      %   state_vector: A matrix of dimension 'n x number_of_states'
      % OPTIONAL:
      %   y_labels    : Cell array containing y-axis labels

      % Create figure or overwrite it
      obj.state_plt = obj.set_fig(obj.state_plt);

      % Set figure properties
      set(obj.state_plt, 'name', 'State Plot', 'color', 'w', ...
          'units', 'normalized');

      if ~strcmp(get(0, 'DefaultFigureWindowStyle'), 'docked')
        set(obj.state_plt, 'outerposition', [0.5 0 0.5 1]);
      end

      % Number of subplots
      plt_num = size(state_vector, 2);

      % Figure heights and top placement point
      height = round(1/plt_num, 3);
      old_top = 1-height;

      % Matrix containing color order
      colors = get(gca, 'ColorOrder');

      j = 1;
      for i = 1:plt_num
        % Create subplot and plot
        subplot(plt_num, 1, j)
        plot(t, state_vector(:, i), 'color', colors(i, :), 'LineWidth', 1.5);
        grid on

        % Apply y-labels if given
        if exist('y_labels', 'var')
          ylabel(['$' y_labels{i} '$'], 'interpreter', 'latex', 'FontSize', 18)
        end

        % Remove xticks for all except the last subplot
        if i < plt_num
          xticklabels([]);
        end

        % Maximize plots (in a very hacky way...)
        p    = get(gca, 'Position');
        p(1) = 0.06; % Distance from the left border
        p(3) = 0.92; % Distance from the right border
        set(gca, 'Position', p);

        op      = get(gca, 'OuterPosition');
        op(2)   = old_top+i*0.009; % Top placement
        op(4)   = 0.95/plt_num;    % Height of plot
        old_top = old_top - height;
        set(gca, 'OuterPosition', op);

        j = j+1;
      end
    end

  end % public methods


  methods (Access = private)

    function debugplt1_internal_f(obj, filepath)
      % Reads file and calls 'debugplt1_internal'.

      % Format of the columns
      formatSpec = '%f %f %f %f %f %f %f %f';
      %
      fileID = fopen(filepath, 'r');

      dataArray = textscan(fileID, formatSpec, 'HeaderLines', 1);

      % Assigning data
      t     = dataArray{1};
      rot_x = dataArray{2};
      rot_y = dataArray{3};
      sta_x = dataArray{4};
      sta_y = dataArray{5};
      theta = dataArray{6};
      fn    = dataArray{7};
      d     = dataArray{8};

      % Call function
      obj.debugplt1_internal(t, rot_x, rot_y, sta_x, sta_y, theta, fn, d)
    end


    function debugplt1_internal(obj, t, rot_x, rot_y, sta_x, sta_y, theta, ...
                                fn, d)
      % Creates a debug plot.

      % Create/append a debug object
      obj.dbugplots{size(obj.dbugplots, 2) + 1} = ...
        Debug2(t, rot_x, rot_y, sta_x, sta_y, theta, fn, d);

    end

  end % private methods


  methods (Static, Access = private)

    function fig_obj = set_fig(fig_handle)
      % Creates a figure or focusses it and thereby overwrites it
      if isempty(fig_handle) || ~isvalid(fig_handle)
        fig_obj = figure();
      else
        set(0, 'CurrentFigure', fig_handle)
        clf
        fig_obj = fig_handle;
      end

    end

  end % static/private methods
end % class
