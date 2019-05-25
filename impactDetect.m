function [value, isterminal, direction] = impactDetect(~, y, s, dir)
% 'impactDetect' used by the ode solvers to halt integration upon rotor-stator
% impact.
%
% INPUTS:
%   y   : State vector
%   s   : Rotor system object
%   dir : Direction of the rotor d = 1 if outgoing and d = -1 if ingoing
%

  value      = s.calc_indent(y);
  isterminal = 1;
  direction  = dir;
end

