classdef Contactmodels < handle
% 'Contactmodel' is a class containing all properties relevant to contact
% models between a rotor and a stator.
% The resulting object is pass-by-reference.

  properties (SetAccess = protected)
    name % Name of the contact model
  end

  properties
    mu_k            % Friction coefficient
  end

  methods
    function obj = Contactmodels()
    % Constructor function.
    % INPUT:
    %   xx: ...
    end

  end

  methods (Access = protected)
    function print_name(obj)
    % 'print_name' displays the current contact model in the console.
      fprintf('Using the %s model.\n', obj.name)
    end
  end
end
