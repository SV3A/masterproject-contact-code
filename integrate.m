function [t_total,y_total,s_total,te_total,crossings] = ...
         integrate(tspan, y_0, s, cmod, app)
  
  % Solver options
  options_ode45 = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, ...
                         'Events', @(t,y) impact_detect(t, y, s, 1));
  options_ode15 = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, ...
                         'Events', @(t,y) impact_detect(t, y, s, -1));

  t_total = 0;
  s_total = 0;
  y_total = [];
  te_total = [];
  crossings = -1;
  loc_tst = tspan(1);

  while t_total(end) ~= tspan(2)
    crossings = crossings + 1;
    
    gap = s.calc_gap(y_0);
    
    if gap < 0
      contact_state = 0;
      [t,y,te,ye,ie] =  ode45(@(t,y) dydt(t,y,s,cmod, contact_state), ...
                        [loc_tst,tspan(2)], y_0, options_ode45);
    else
      contact_state = 1;
      [t,y,te,ye,ie] = ode15s(@(t,y) dydt(t,y,s,cmod,contact_state), ...
                       [loc_tst,tspan(2)], y_0, options_ode15);
    end
    
    % Collect results
    t_total = [ t_total(1:end-1)  ; t ];
    y_total = [ y_total(1:end-1,:); y ];
    s_total = [ s_total(1:end-1)  ; contact_state*ones(length(t),1) ];
    te_total = [te_total; te];
    
    % Assign new initial conditions
    loc_tst = t(end);
    y_0 = y(end,:);
    t_total(end)
    
    % Check if run from GUI
    if nargin > 4
      if ~app.isrunning
        return
      end
      
      app.status_line_lable.Text = ['Running...  '...
        num2str(round((t_total(end)/tspan(2))*100),2),'% completed.'];
      drawnow
    end
    
  end
end
