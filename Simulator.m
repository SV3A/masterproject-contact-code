classdef Simulator < handle
  properties (Access = public)
    % Solution properties
    time = 0           % time
    solution = []      % solution
    event_times = []   % event times
    contact_states = 0 % contact state vector
    y_0 = zeros(14,1); % initial condition vector

    % Model objects
    s                  % rotor-stator system
    cmod               % contact model

    % Model properties
    xi
    m0
    e
    fric_mod

    % Derived parameters
    r_OC
    s_OC
    F_c
    fn
    d
  end

  methods

    function obj = Simulator()
    % Constructor

    end

    function solve(obj, tspan)
      % Initiate system object with parameters 'xi', 'm0' and 'e'
      obj.s = Rotorsystem(obj.xi, obj.m0, obj.e);

      % Define the contact model
      obj.cmod = Nikravesh(obj.fric_mod, obj.s.r_s, obj.s.r_r);

      % Solver options
      options_ode45 = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', 1e-3,...
                             'Events', @(t,y) impactDetect(t, y, obj.s, 1));
      options_ode15 = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, ...
                             'Events', @(t,y) impactDetect(t, y, obj.s, -1));

      loc_tst = tspan(1); % Integration time starting point

      tic
      while obj.time(end) ~= tspan(2)

        indent = obj.s.calc_indent(obj.y_0);

        if indent <= 0
          contact_state = 0;
          [t, y, te] =  ode45(@(t,y) dydt(t, y, obj.s, obj.cmod, ...
            contact_state), [loc_tst, tspan(2)], obj.y_0, options_ode45);
        else
          contact_state = 1;
          [t, y, te] = ode15s(@(t,y) dydt(t, y, obj.s, obj.cmod, ...
            contact_state), [loc_tst, tspan(2)], obj.y_0, options_ode15);
        end

        % Collect results
        obj.time  = [ obj.time(1:end-1)  ; t ];
        obj.solution  = [ obj.solution(1:end-1,:); y ];
        obj.contact_states  = [ obj.contact_states(1:end-1); ...
          contact_state*ones(length(t),1) ];
        obj.event_times = [obj.event_times; te];

        % Assign new initial conditions
        loc_tst = t(end);
        obj.y_0 = y(end,:);

        obj.time(end)
      end
      toc

      fprintf('%i perimeter crossings detected\n', length(obj.event_times))
    end


    function postprocess(obj)
    % 'postprocess' calculate the forces associated with a given solution.

      % Check if solution
      if isempty(obj.solution)
        error('No solution present.')
      end

      % Init
      obj.r_OC = zeros( 3, length(obj.time) );
      F_cxs    = zeros( length(obj.time), 1 );
      F_cys    = zeros( length(obj.time), 1 );
      obj.d        = zeros( length(obj.time), 1 );

      % Build orbit in I, get contact forces, retrieve the contact angle and get
      % relative radial velocity
      for i = 1:length(obj.time)
        y_i = obj.solution(i,:)';

        T_gamma = obj.s.T_gam(obj.solution(i, 1));
        T_beta  = obj.s.T_bet(obj.solution(i, 3));
        obj.r_OC(:,i) = T_gamma' * (T_beta'*[0; 0; obj.s.l_OC]);

        if obj.contact_states(i) == 0, state = 0; else, state = 1;  end

        [F_cxs(i), F_cys(i)] = contactForce(y_i, obj.s, obj.cmod, state);

        obj.d(i) = obj.s.calc_indent(y_i);

      end

      obj.F_c = [F_cxs'; F_cys'; zeros(1, size(obj.solution, 1))];

      % Calculate resultant radial contact force
      obj.fn = sqrt(F_cxs.^2 + F_cys.^2);

      % Stator centre in the plane of contact
      obj.s_OC = [
        obj.solution(:, 7)'
        obj.solution(:, 9)'
        zeros(size(obj.solution, 1)) ];
    end

    function export(obj, export_type)
    % 'export' handle for function to export the solution to a text file.

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
