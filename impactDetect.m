function [value, isterminal, direction] = impactDetect(~, y, s, dir, term)
% 'impactDetect' used by the ode solvers to halt integration upon rotor-stator
% impact.
%
% INPUTS:
%   y   : State vector
%   s   : Rotor system object
%   dir : Direction of the rotor d = 1 if outgoing and d = -1 if ingoing
%   term: Terminate on trigger 0|1
%

  value      = s.calc_indent(y);
  isterminal = term;
  direction  = dir;
end

