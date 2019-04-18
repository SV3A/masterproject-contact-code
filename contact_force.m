function [F_cx, F_cy, delta, delta_d] = contact_force(y, s, cmod, state)
% 'contact_force' computes the contact force and related parameters.
%
% INPUT:
%   y   : State vector, the substitution is given as:
%         y = [gamma, gamma_d, beta, beta_d, theta, theta_d,
%              x_ih, x_ih_d, y_ih, y_ih_d, x_mh, x_mh_d, y_mh, y_mh_d]^T
%     s   : Rotor system object
%     cmod: The contact model object
%    state: Direction, "1" means outgoing while "-1" is ingoing
% OUTPUT:
%     F_c*: Force component
%    delta: Indentation
%  delta_d: Indentation rate
%

  if s.calc_gap(y) < 0
  % If the "gap" is negative the contact force and indentation is zero, while
  % the initial relative impact velocity should be calculated
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
    delta = s.calc_gap(y);

    % Penetration rate
    delta_d = s.pen_rate(y);

    % Contact angle
    alpha = s.contact_ang(y);

    % Normal force
    Fn = cmod.calc_fn(delta, delta_d);

    % Contact forces in the inertial coordinate system
    F = Fn*[cos(alpha) - cmod.mu_k*sin(alpha)*sign(y(6));
            sin(alpha) + cmod.mu_k*cos(alpha)*sign(y(6));
            0];

    F_cx = F(1);
    F_cy = F(2);

    delta_d = cmod.delta_d_init;
  end
end

