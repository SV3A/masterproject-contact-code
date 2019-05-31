classdef Plottools < handle
  % 'Plottools' is a class which defines a common plotting module.

  properties (Access = private)
    dbugplots = {};  % Cell array containing debug-plot objects
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
end % class
