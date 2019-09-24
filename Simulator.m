classdef Simulator < handle
  % 'Simulator' is a class defining a interface for simulating the rotor-stator
  % system.

  properties
    % Solution properties

    time           % Time
    solution       % Solution
    event_times    % Event times
    contact_states % Contact state vector
    y_0            % Initial condition vector

    % Model objects

    s         % Rotor-stator system
    cmod      % Contact model

    % Model properties

    fric_mod  % Friction model
    clearance % Static clearance between the rotor and stator

    % External excitation magnet

    mag_enabled; % Boolean for knowing if the magnet has been enabled
    mag_app_t; % Absolute time, rel. to sim. start, of when to apply the magnet
    mag_app_angle; % Angle to sync the magnet application with
    mag_forcedata; % Force data nx2 vector containing time [s] force [N]

    % Derived parameters

    r_OC      % Position vector of the rotor centre in the contact plane in I
    r_OD      % Position vector of the rotor centre in the disc plane in I
    s_OC      % Position vector of the stator centre in the contact plane in I
    F_c       % Contact force vector
    fn        % Magnitude of the normal force
    d         % Radial indentation

    % Solver tolerances

    o45_reltol % Relative tolerance for the ode45 solver
    o45_abstol % Absolute tolerance for the ode45 solver
    o15_reltol % Relative tolerance for the ode15s solver
    o15_abstol % Absolute tolerance for the ode15s solver
  end


  methods

    function obj = Simulator()
      % Constructor function.

      % Set default solver tolerances
      obj.o45_reltol = 1e-9;
      obj.o45_abstol = 1e-9;
      obj.o15_reltol = 1e-9;
      obj.o15_abstol = 1e-9;

      % Initial condition
      obj.y_0 = zeros(14, 1);

      obj.mag_enabled = false;
    end


    function set_magnet(obj, varargin)
      % Enables the external excitation magnet and sets its properties.
      %
      % INPUT:
      %   time : the rough point in time when to apply the magnet force
      %   angle: (optional) the angle at which to apply the magnet
      %   force: a nx2 vector containing the local time and force of the magnet

      if ~obj.mag_enabled

        % Handle optional arguments (note nargin is +1 due to object)
        if nargin < 3
          error('ERROR Not enough inputs supplied')
        elseif nargin > 3
          obj.mag_app_angle = varargin{2};
          obj.mag_forcedata = varargin{3};
        else
          obj.mag_app_angle = 0;
          obj.mag_forcedata = varargin{2};
        end

        obj.mag_app_t = varargin{1};

        if size(obj.mag_forcedata, 2) ~= 2
          error('ERROR Magnet force data not a vector of dimension nx2')
        end

        obj.mag_enabled = true;
      else
        warning("Magnet already enabled, I'm ignoring this.")
      end
    end


    function ssolve(obj, tspan)
      % Performs the time integration using only the 15s solver.

      % Init/clear
      obj.time           = 0;
      obj.solution       = [];
      obj.event_times    = [];
      obj.contact_states = 0;

      % Initiate system object
      if obj.mag_enabled
        obj.s  = Rotorsystem(obj.mag_app_t, obj.mag_app_angle, ...
                             obj.mag_forcedata);
      else
        obj.s  = Rotorsystem();
      end
      obj.cmod = Nikravesh(obj.fric_mod, obj.s.r_s, obj.s.r_r);

      % Solver options
      options_ode15 = odeset('RelTol', obj.o15_reltol, ...
                             'AbsTol', obj.o15_abstol, ...
                             'Events', @(t,y) impactDetect(t, y, obj.s, 0, 0));

      tic
      [t, y, te] = ode15s(@(t,y) dydt(t, y, obj.s, obj.cmod, 0), tspan, ...
                          obj.y_0, options_ode15);
      toc

      % Collect results
      obj.time     = [obj.time(1:end-1)      ; t];
      obj.solution = [obj.solution(1:end-1,:); y];
      obj.event_times = te;

      % Find contact states
      for i = 1:size(y, 1)
        if obj.s.calc_indent(y(i,:)) <= 0
          obj.contact_states(i) = 0;
        else
          obj.contact_states = 1;
        end
      end

      fprintf('%i perimeter crossings detected\n', length(obj.event_times))
    end

    function solve(obj, tspan)
      % Performs the time integration.

      % Init/clear
      obj.time           = 0;
      obj.solution       = [];
      obj.event_times    = [];
      obj.contact_states = 0;

      % Initiate system objects
      if obj.mag_enabled
        obj.s  = Rotorsystem(obj.mag_app_t, obj.mag_app_angle, ...
                             obj.mag_forcedata);
      else
        obj.s  = Rotorsystem();
      end
      obj.cmod = Nikravesh(obj.fric_mod, obj.s.r_s, obj.s.r_r);

      % Solver options
      options_ode45 = odeset('RelTol', obj.o45_reltol, ...
                             'AbsTol', obj.o45_abstol, 'MaxStep', 1e-3,...
                             'Events', @(t,y) impactDetect(t, y, obj.s, 1, 1));

      options_ode15 = odeset('RelTol', obj.o15_reltol, ...
                             'AbsTol', obj.o15_abstol, ...
                             'Events', @(t,y) impactDetect(t, y, obj.s, -1, 1));

      loc_tst = tspan(1); % Integration time starting point
      y0 = obj.y_0;       % Initial condition

      tic
      while obj.time(end) ~= tspan(2)

        indent = obj.s.calc_indent(y0);

        if indent < 0 || ( indent == 0 && obj.s.pen_rate(y0) < 0 )
          contact_state = 0;
          [t, y, te] =  ode45(@(t,y) dydt(t, y, obj.s, obj.cmod, ...
            contact_state), [loc_tst, tspan(2)], y0, options_ode45);
        else
          contact_state = 1;
          [t, y, te] = ode15s(@(t,y) dydt(t, y, obj.s, obj.cmod, ...
            contact_state), [loc_tst, tspan(2)], y0, options_ode15);
        end

        % Collect results
        obj.time     = [obj.time(1:end-1)      ; t];
        obj.solution = [obj.solution(1:end-1,:); y];
        obj.contact_states  = [obj.contact_states(1:end-1); ...
                               contact_state*ones(length(t),1)];
        obj.event_times = [obj.event_times; te];

        % Assign new initial conditions
        loc_tst = t(end);
        y0      = y(end,:);

        % Print feedback
        fprintf('t_n = %f s\n', obj.time(end));
      end
      toc

      fprintf('%i perimeter crossings detected\n', length(obj.event_times))
    end


    function postprocess(obj)
      % Calculates the forces associated with a given solution.

      % Check if solution is present
      if isempty(obj.solution)
        error('No solution present.')
      end

      % Init/clear
      F_cxs    = zeros( length(obj.time), 1 );
      F_cys    = zeros( length(obj.time), 1 );
      obj.r_OC = zeros( 3, length(obj.time) );
      obj.r_OD = zeros( 3, length(obj.time) );
      obj.s_OC = zeros( 3, length(obj.time) );
      obj.F_c  = zeros( 3, length(obj.time) );
      obj.fn   = zeros( length(obj.time), 1 );
      obj.d    = zeros( length(obj.time), 1 );

      % Build orbit in I, get contact forces, retrieve the contact angle and get
      % relative radial velocity
      for i = 1:length(obj.time)
        y_i = obj.solution(i,:)';

        T_gamma = obj.s.T_gam(obj.solution(i, 1));
        T_beta  = obj.s.T_bet(obj.solution(i, 3));
        obj.r_OC(:,i) = T_gamma' * (T_beta'*[0; 0; obj.s.l_OC]);
        obj.r_OD(:,i) = T_gamma' * (T_beta'*[0; 0; obj.s.l_OD]);

        if obj.contact_states(i) == 0, state = 0; else, state = 1;  end

        [F_cxs(i), F_cys(i)] = contactForce(y_i, obj.s, obj.cmod, state);

        obj.d(i) = obj.s.calc_indent(y_i);
      end

      obj.F_c = [F_cxs'; F_cys'; zeros(1, size(obj.solution, 1))];

      % Calculate resultant radial contact force
      obj.fn = sqrt(F_cxs.^2 + F_cys.^2);

      % Stator centre in the plane of contact
      obj.s_OC = [obj.solution(:, 7)';
                  obj.solution(:, 9)';
                  zeros(1, size(obj.solution, 1))];

      obj.clearance = obj.s.cl;
    end


    function export(obj, export_type)
      % A handle for function to export the solution to a text file.

      % Define parameter list
      if strcmp(export_type, 'basic')
        par_list = {'t', 'rotor_x', 'rotor_y',  'stator_x', 'stator_y', ...
                    'theta', 'Fn', 'delta'};

        value_vector = [obj.time'; obj.r_OC(1, :); obj.r_OC(2, :); ...
                        obj.s_OC(1, :); obj.s_OC(2, :); obj.solution(:, 5)'; ...
                        obj.fn'; obj.d'];

      elseif strcmp(export_type, 'all')
        %par_list = {'t', 'rotor_x', 'rotor_y',  'stator_x', 'stator_y', ...
                    %'theta', 'Fn', 'delta', 'x', 'y', 'z'};

        %value_vector = [value_vector; x; y; z;];
      else
        fclose(fileID); error('Unknown export type.')
      end

      % Call external function
      export_values(par_list, value_vector)
    end

  end % methods
end % class
