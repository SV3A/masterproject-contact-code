classdef Contactmodels < handle
% 'Contactmodel' is a class containing all properties relevant to contact
% models between a rotor and a stator.
% The resulting object is pass-by-reference.

  properties (SetAccess = protected)
    name % Name of the contact model
  end

  properties
    mu_k            % Friction coefficient
    v_0  = 0.0001;  % [m/s]  lower velocity tolerance
    v_1  = 0.0002;  % [m/s]  upper velocity tolerance
    calc_ff
  end

  methods
    function obj = Contactmodels(friction_model)
    % Constructor function.
    % INPUT:
    %   friction_model: Parameter defining which friction model to use

      obj.mu_k = 0.2;

      % Assign function handle for the friction model
      if strcmp(friction_model, 'ambrosio')
        obj.calc_ff = @obj.calc_ff_ambrosio;
      elseif strcmp(friction_model, 'none')
        obj.calc_ff = @(~, ~) 0;
      else
        error('Invalid friction model')
      end
    end

    function Ff = calc_ff_ambrosio(obj, Fn, vt_rel)
    % 'calc_ff_ambrosio' implements the friction model set forth in the
    % article "Influence of the contact–impact force model on the dynamic
    % response of multi-body systems"
    % In: Proc. IMechE, Year 2006, Volume 220.
    % By: Flores, P. and Ambrósio, J. and Claro, J. C. P. and Lankarani, H. M.
    %
    % INPUT:
    %   Fn: Normal force (scalar)
    %   vt_rel: Relative tangential velocity

      % Dynamical correction coefficient
      if (obj.v_1 < norm(vt_rel))
        n_d = 1;
      elseif obj.v_0 < norm(vt_rel) & norm(vt_rel) <= obj.v_1
        n_d = (norm(vt_rel) - obj.v_0)/(obj.v_1 - obj.v_0);
      elseif norm(vt_rel) <= obj.v_0
        n_d = 0;
      end

      Ff = Fn * obj.mu_k * sign(vt_rel) * n_d;
    end
  end

  methods (Access = protected)
    function print_name(obj)
    % 'print_name' displays the current contact model in the console.
      fprintf('Using the %s model.\n', obj.name)
    end
  end
end
