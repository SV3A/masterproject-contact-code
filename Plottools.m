classdef Plottools < handle
  % 'Plottools' is a class which defines a common plotting module.

  properties (Access = private)
    dbugplots = {};  % Cell array containing debug-plot objects

    orbit_plt % Orbit plot figure handle
  end


  properties (Access = public)
    %
  end


  methods (Access = public)

    function obj = Plottools()
      % Constructor function
    end


    function debugplot(obj, varargin)
      % 'debugplot' is a public function handle for creating a debug plot.  The
      % function takes 1 input (filename) or 8 (see 'debugplt1_internal').

      % If only one arg is given it should be a file to read inputs from
      if length(varargin) == 1
        obj.debugplt1_internal_f(varargin{1})
      elseif length(varargin) == 8
        obj.debugplt1_internal(varargin{:})
      else
        error('Too many arguments given')
      end
    end


    function orbit(obj, rot_x, rot_y, clearance)
      % plots the rotor orbit inside the stator boundary (clearance circle).

      % Create figure or overwrite it
      obj.orbit_plt = obj.set_fig(obj.orbit_plt);

      % Set figure properties
      set(obj.orbit_plt, 'name', 'Orbit Plot', ...
          'units', 'normalized', 'outerposition', [0.5 0 0.5 1]);

      % Rotor orbit (converted to [mm])
      plot(rot_x * 1e3, rot_y * 1e3, 'b'); grid on; hold on

      % Clearance circle
      plot(clearance * cos(linspace(0, 2*pi)) * 1e3, ...
           clearance * sin(linspace(0, 2*pi)) * 1e3, 'r')
      xlabel('Orbit [mm]')
      axis equal
    end

  end % public methods


  methods (Access = private)

    function debugplt1_internal_f(obj, filepath)
      % reads file and calls 'debugplt1_internal'.

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
      % 'debugplt1_internal' creates a debug plot.

      % Create/append a debug object
      obj.dbugplots{size(obj.dbugplots, 2) + 1} = ...
        Debug2(t, rot_x, rot_y, sta_x, sta_y, theta, fn, d);

    end

  end % private methods


  methods (Static, Access = private)

    function fig_obj = set_fig(fig_handle)
      % creates a figure or focusses it and thereby overwrites it
      if isempty(fig_handle)
        fig_obj = figure();
      else
        set(0, 'CurrentFigure', fig_handle)
        clf
        fig_obj = fig_handle;
      end

    end

  end % static/private methods
end % class
