function y_dot = dydt(t, y, s, cmod, state)
% 'dydt' serves the equation of motion as a system of six first order equations.
%
% INPUT:
%   y     : State vector
%   s     : Rotor system object
%   cmod  : The contact model object
%   state : Contact indicator (0|1)
% OUTPUT:
%   y_dot: Accelerations and velocities
%

  l_OD   = s.l_OD;
  l_OM   = s.l_OM;
  l_OC   = s.l_OC;
  l_OE   = s.l_OE;

  I_xx   = s.I_xx;
  I_yy   = s.I_yy;
  I_zz   = s.I_zz;

  k_xx   = s.k_xx;
  k_yy   = s.k_yy;
  k_xy   = s.k_xy;
  k_yx   = s.k_yx;
  k_ft1  = s.k_ft1;
  k_ft2  = s.k_ft2;
  k_vg   = s.k_vg;
  k_hg   = s.k_hg;

  m0     = s.m0;
  e_x    = s.e_x;
  e_y    = s.e_y;

  m_ih   = s.m_ih;
  m_mh   = s.m_mh;

  d_vg   = s.d_vg;
  d_hg   = s.d_hg;
  d_yy   = s.d_yy;

  if t < 2
    d_xx = 100;
    d_yy = 100;
  end

  % Damping switch to account for anisotropic damping at standstill
  if abs(y(6)) > 0
    d_xx  = s.d_xx;
  else
    d_xx  = s.d0_xx;
  end

  % Calculate contact forces
  [F_cx, F_cy] = contactForce(y, s, cmod, state);

  % Evaluate external magnet force if enabled
  if s.mag_enabled
    [F_excx, F_excy] = s.magnetForce(t, y);
  else
    F_excx = 0;
    F_excy = 0;
  end


% Don't edit anything after this as it may be overwritten by Maple!
% * * * * *
  y_dot(1)  = y(2);
  y_dot(2)  = 0.1e1 / I_xx / cos(y(3)) * (-cos(y(5)) * e_y * l_OD * m0 * y(6) ^ 2 - sin(y(5)) * e_x * l_OD * m0 * y(6) ^ 2 - y(1) * cos(y(1)) * k_yy * l_OM ^ 2 + cos(y(1)) * y(3) * k_yx * l_OM ^ 2 - cos(y(1)) * d_yy * l_OM ^ 2 * y(2) + I_xx * sin(y(3)) * y(4) * y(2) + sin(y(3)) * y(2) * I_yy * y(4) - I_zz * sin(y(3)) * y(4) * y(2) + l_OC * cos(y(1)) * F_cy - l_OE * cos(y(1)) * F_excy - I_zz * y(4) * y(6));
  y_dot(3)  = y(4);
  y_dot(4)  = -(-y(1) * sin(y(3)) * sin(y(1)) * k_yy * l_OM ^ 2 + sin(y(3)) * sin(y(1)) * y(3) * k_yx * l_OM ^ 2 - sin(y(3)) * sin(y(1)) * d_yy * l_OM ^ 2 * y(2) - cos(y(5)) * e_x * l_OD * m0 * y(6) ^ 2 + sin(y(5)) * e_y * l_OD * m0 * y(6) ^ 2 - y(1) * cos(y(3)) * k_xy * l_OM ^ 2 + sin(y(3)) * y(2) ^ 2 * I_xx * cos(y(3)) - I_zz * cos(y(3)) * sin(y(3)) * y(2) ^ 2 + cos(y(3)) * y(3) * k_xx * l_OM ^ 2 + cos(y(3)) * y(4) * d_xx * l_OM ^ 2 + F_cy * sin(y(3)) * sin(y(1)) * l_OC - F_excy * sin(y(3)) * sin(y(1)) * l_OE - I_zz * cos(y(3)) * y(2) * y(6) + F_cx * cos(y(3)) * l_OC - F_excx * cos(y(3)) * l_OE) / I_yy;
  y_dot(5)  = y(6);
  y_dot(6)  = -(-I_zz * sin(y(3)) * cos(y(5)) * e_y * l_OD * m0 * y(6) ^ 2 - I_zz * sin(y(3)) * sin(y(5)) * e_x * l_OD * m0 * y(6) ^ 2 - y(1) * I_zz * sin(y(3)) * cos(y(1)) * k_yy * l_OM ^ 2 + I_zz * sin(y(3)) * cos(y(1)) * y(3) * k_yx * l_OM ^ 2 - I_zz * sin(y(3)) * cos(y(1)) * d_yy * l_OM ^ 2 * y(2) - cos(y(3)) ^ 2 * y(2) * y(4) * I_xx ^ 2 + cos(y(3)) ^ 2 * y(2) * y(4) * I_xx * I_yy + cos(y(3)) ^ 2 * y(2) * y(4) * I_xx * I_zz + I_xx * I_zz * sin(y(3)) ^ 2 * y(4) * y(2) + I_yy * I_zz * sin(y(3)) ^ 2 * y(4) * y(2) - I_zz ^ 2 * sin(y(3)) ^ 2 * y(4) * y(2) + F_cy * I_zz * sin(y(3)) * cos(y(1)) * l_OC - F_excy * I_zz * sin(y(3)) * cos(y(1)) * l_OE - I_zz ^ 2 * sin(y(3)) * y(4) * y(6)) / I_xx / cos(y(3)) / I_zz;
  y_dot(7)  = y(8);
  y_dot(8)  = (-d_vg * y(8) + d_vg * y(12) - k_vg * y(7) + k_vg * y(11) + F_cx) / m_ih;
  y_dot(9)  = y(10);
  y_dot(10) = (-2 * k_ft1 * y(9) + 2 * k_ft1 * y(13) + F_cy) / m_ih;
  y_dot(11) = y(12);
  y_dot(12) = (d_vg * y(8) - d_vg * y(12) - 2 * k_ft2 * y(11) + k_vg * y(7) - k_vg * y(11)) / m_mh;
  y_dot(13) = y(14);
  y_dot(14) = -1 / m_mh * (d_hg * y(14) - 2 * k_ft1 * y(9) + 2 * k_ft1 * y(13) + k_hg * y(13));

  y_dot = y_dot';

end
% Maple part last updated: 22-09-2019 15:31
