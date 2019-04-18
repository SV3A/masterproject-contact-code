function [value, isterminal, direction] = impact_detect(~, y, s, dir)
% 'impact detect' used by the ode solvers to halt integration upon rotor-stator
% impact.
%
% INPUTS:
%   y  : State vector, the substitution is given as:
%        y = [gamma, gamma_d, beta, beta_d, theta, theta_d, 
%             x_ih, x_ih_d, y_ih, y_ih_d, x_mh, x_mh_d, y_mh, y_mh_d]^T
%   s  : Rotor system object
%   dir: Direction of the rotor d = 1 if outgoing and d = -1 if ingoing
%

  value = s.calc_gap(y);
  isterminal = 1;
  direction  = dir;
end
