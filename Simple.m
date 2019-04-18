classdef Simple < Contactmodels
% 'Simple' is a subclass of 'Contactmodels' and defines the contact model
% only dependent on penetration.

  properties (Constant)
    k = 1e6   % Stator material stiffness
  end
  
  methods
    function obj = Simple()
    % Constructor function.
      obj.name = "Simple";
      obj.print_name
    end
    
    function Fn = calc_fn(obj, d, ~)
    % 'calc_fn' calculates the normal force
      Fn = obj.k*d;
    end
  end
end
