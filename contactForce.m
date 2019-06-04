function [F_cx, F_cy, delta, delta_d] = contactForce(y, s, cmod, state)
% 'contactForce' computes the contact force and related parameters.
%
% INPUT:
%   y     : State vector
%   s     : Rotor system object
%   cmod  : The contact model object
%   state : Contact indicator (0|1)
% OUTPUT:
%   F_c*    : Force components
%   delta   : Indentation
%   delta_d : Indentation rate
%

  if s.calc_indent(y) <= 0
  % If the "indent" is negative the contact force and indentation is zero, while
  % the initial relative impact velocity should be calculated

    % Check which solver is calling the parent function, this is important since
    % delta_d_init cannot change, when the event function is only "testing" the
    % limit
    if state == 0
      cmod.delta_d_init = s.pen_rate(y);
      delta_d = cmod.delta_d_init;
    else
      delta_d = 0;
    end
    F_cx = 0;
    F_cy = 0;
    delta = 0;

  else
    % Penetration
    delta = s.calc_indent(y);

    % Penetration rate
    delta_d = s.pen_rate(y);

    % Relative velocity between the rotor- and stator surface (that is the
    % tangential component)
    vt_rel = s.tan_rel_velocity(y);

    % Normal- and friction force
    Fn = cmod.calc_fn(delta, delta_d);
    Ff = cmod.calc_ff(Fn, vt_rel);

    % Contact angle
    alpha = s.contact_ang(y);

    % Contact forces in the inertial coordinate system
    F = s.T_the(alpha)' * [Fn; Ff; 0];

    F_cx = F(1);
    F_cy = F(2);

    % delta_d is overwritten here for debug purposes in 'debug1'
    delta_d = cmod.delta_d_init;
  end
end

