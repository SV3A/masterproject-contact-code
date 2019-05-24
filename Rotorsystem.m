classdef Rotorsystem < handle
% 'Rotorsystem' is a class containing all properties relevant to a particular
% rotor-dynamical system.  The resulting object is pass-by-reference.
% The state vector substitution is given as:
%   y = [gamma, gamma_d, beta, beta_d, theta, theta_d,
%        x_ih, x_ih_d, y_ih, y_ih_d, x_mh, x_mh_d, y_mh, y_mh_d]^T

  properties (Constant)
    % Gravity constant [m/s^2]
    g = 9.81;

    r_r = 25e-3/2;    % Rotor radius [m]
    r_s = 29e-3/2;    % Stator radius [m]

    l_OM = 171.7e-3;  % Position vector, pivot point to the PMB [m]
    l_OG = 195.9e-3;  % Position vector, pivot point to centre of gravity [m]
    l_OD = 259.7e-3;  % Position vector, pivot point to the disc [m]
    l_OC = 413.7e-3;  % Position vector, pivot point to contact point [m]

    m_tot = 1.486;    % Shaft mass (including discs) [kg]

    I_xx = 73399360 * 1e-9; % Mass moment of inertia component [kg*m^2]
    I_yy = 73399360 * 1e-9; % Mass moment of inertia component [kg*m^2]
    I_zz = 981365 * 1e-9;   % Mass moment of inertia component [kg*m^2]

    K_mx = 3.4675*10^4;     % Stiffness of the magnetic bearing in x [N/m]
    K_my = 3.4675*10^4;     % Stiffness of the magnetic bearing in y [N/m]

    l_OIH = 413.7e-3;   % Position vector, pivot point to inner house [m]
    l_OMH = 413.7e-3;   % Position vector, pivot point to middle [m]

    m_ih = 1.7;     % Inner house mass [kg]
    m_mh = 7.17;    % Middle house mass [kg]

    K_ft = 8.88e7;  % Stiffness of the force transducer [N/m]

    K_vg = 1.23e7;  % Stiffness of vertical beams [N/m]
    K_hg = 2.40e7;  % Stiffness of horizontal beams [N/m]

    D_vg = 238;     % Damping coefficient of vertical beams [Ns/m]
    D_hg = 1576;    % Damping coefficient of horizontal beams [Ns/m]

    x0 = 0;   % Initial housing offset in x [m]
    y0 = 0;   % Initial housing offset in y [m]

  end

  properties
    m0      % Unbalance mass
    e       % Eccentricity
    xi      % Damping ratio [-]
    D_x     % Damping coefficient in x [N*s/m]
    D_y     % Damping coefficient in y [N*s/m]
    cl      % Clearance [m]
  end

  methods
    function obj = Rotorsystem(xi, m0, e)
    % Constructor function.
    % INPUT:
    %   xi: damping ratio
    %   m0: unbalance mass
    %   e : unbalance eccentricity

      obj.m0 = m0;
      obj.e  = e;
      obj.xi = xi;

      % Calculate damping factor
      obj.D_x = 2*obj.xi*sqrt( obj.K_mx*obj.m_tot );
      obj.D_y = 2*obj.xi*sqrt( obj.K_my*obj.m_tot );

      % Calculate clearance
      obj.cl = obj.r_s - obj.r_r;
    end

    function T_gamma = T_gam(~, gamma)
    % Transformation matrix from I  to B1
      T_gamma = [1       0          0
                 0  cos(gamma) sin(gamma)
                 0 -sin(gamma) cos(gamma)];
    end

    function T_beta = T_bet(~, beta)
    % Transformation matrix from B1 to B2
      T_beta = [cos(beta) 0 -sin(beta)
                    0     1      0
                sin(beta) 0  cos(beta)];
    end

    function T_theta = T_the(~, theta)
    % Transformation matrix from B2 to B3
      T_theta = [ cos(theta) sin(theta) 0
                 -sin(theta) cos(theta) 0
                       0          0     1];
    end

    function r_OR = rot_centrepos(obj, y)
    % 'r_OR(y)' gets position of the rotor centre in the contact plane in the
    % inertial coordinate system.
    %
    % INPUT:
    %   y: State vector.

      r_OR = obj.T_gam(y(1))' * (obj.T_bet(y(3))' * [0; 0; obj.l_OC]);
    end

    function gap = calc_gap(obj, y)
    % 'calc_gap(y)' calculates the size of the gap between the rotor and stator.
    %
    % INPUT:
    %   y: State vector.
    %
    % OUTPUT:
    %   gap: The size of the 'gap', where the gap is represented as the rotor
    %   centre position minus the clearance, thus the output is negative when
    %   the rotor is not in contact with the stator, and positive otherwise,
    %   during which the output signifies the penetration of the rotor into the
    %   stator material.
    %
      % Get rotor centre position
      r_OR = obj.rot_centrepos(y);

      % Calculate the effective clearance (or penetration)
      gap = sqrt( (r_OR(1)-y(7))^2 + (r_OR(2)-y(9))^2 ) - obj.cl;
    end

    function alpha = contact_ang(obj, y)
    % 'contact_ang(y)' calculates the contact angle between houses and rotor.
    % INPUT:
    %   y: State vector.

      % Get rotor centre position
      r_OR = obj.rot_centrepos(y);

      alpha = atan2(r_OR(2)-y(9), r_OR(1)-y(7));
    end

    function v_rc = rotor_linvel(obj, y)
    % 'rotor_linvel(y)' calculates the linear velocity of the rotor centre in
    % the inertial coordinate system.
    % INPUT:
    %   y: State vector.

      % Absolute angular velocity of the reference frame B2
      Omega_B2 = obj.T_bet(y(3)) * obj.T_gam(y(1)) * [y(2); 0; 0] + ...
                 obj.T_bet(y(3)) * [0; y(4); 0];

      % Absolute linear velocity of the rotor centre in the contact plane
      v_rc = obj.T_gam(y(1))' * (obj.T_bet(y(3))' * ...
             cross( Omega_B2, [0; 0; obj.l_OC] ));
    end

    function v_rel = abs_rel_velocity(obj, y)
    % 'rel_velocity' calculates the relative velocity between the rotor and
    % the stator in the inertial coordinate system.
    % INPUT:
    %   y: State vector.

      v_rel = obj.rotor_linvel(y) - [y(7); y(9); 0];
    end

    function vt_rel = tan_rel_velocity(obj, y)
    % 'tan_rel_velocity' calculates the tangential compoent of the relative
    % velocity between the rotor and stator.
    % INPUT:
    %   y: State vector.

      v_rel = obj.abs_rel_velocity(y);

      % Get contact angle between rotor and inner house
      alpha = obj.contact_ang(y);

      vt_rel = obj.r_r*y(6) - v_rel(1)*sin(alpha) + v_rel(2)*cos(alpha);
    end

    function [delta_d, v_rel_r] = pen_rate(obj, y)
    % 'pen_rate' calculates the penetration rate.
    % INPUT:
    %   y: State vector.

      % Get contact angle between rotor and inner house
      alpha = obj.contact_ang(y);

      % Unit vector in the radial direction
      e_r = [cos(alpha); sin(alpha); 0];

      % Relative velocity between the rotor and the inner house (stator)
      v_rel = obj.abs_rel_velocity(y);

      % Projection onto the radial direction
      v_rel_r = dot(v_rel, e_r)*e_r;

      % Determine direction of the identation rate
      if sign(v_rel_r(1)) == sign(e_r(1)) && sign(v_rel_r(2)) == sign(e_r(2))
        dir =  1;
      else
        dir = -1;
      end

      % Penetration rate, i.e. the magnitude of the radial velocity
      delta_d = norm( v_rel_r )*dir;
    end

  end
end
