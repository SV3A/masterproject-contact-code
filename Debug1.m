classdef Debug1 < handle
  % 'Debug1' is a class used to create a visual debugging tool.

  properties
    % Solution data
    t_total
    y_total
    te_total
    s_total
    s
    cmod

    % Solution parameters
    pos1
    alphas
    v_rel_r
    F_cxs
    F_cys
    F_r
    deltas
    delta_ds
    delta_d_ini

    % Plot and GUI handles
    mainfig
    rotor
    stator
    rotor_cross
    stator_cross
    rev_mark
    v_vector
    force_plt
    forcein_plt
    info_field
    fplt
    imp_field
    sld_handle
    sens_button
    fwd_button
    bwd_button
  end % properties

  methods (Access = public)
    function obj = Debug1( t_total, y_total, te_total, s_total, s, cmod )
      % Constructor function.
      % INPUT:
      %   t_total  : Time steps
      %   y_total  : Solution step
      %   te_total : Event times
      %   s_total  : Contact states
      %   s        : Rotor system object
      %   cmod     : Contact model object

      obj.t_total  = t_total;
      obj.y_total  = y_total;
      obj.te_total = te_total;
      obj.s_total  = s_total;
      obj.s        = s;
      obj.cmod     = cmod;

      obj.setupSol();
      obj.setupPlots();
      obj.setupGUI();
    end
  end

  methods (Access = private)
    function setupSol(obj)
    % 'setupSol' calculates solution parameters.

      % Init
      obj.pos1     = zeros( 3, length(obj.t_total) );
      obj.alphas   = zeros( length(obj.t_total), 1 );
      obj.v_rel_r  = zeros( 3, length(obj.t_total) );
      obj.F_cxs    = zeros( length(obj.t_total), 1 );
      obj.F_cys    = zeros( length(obj.t_total), 1 );
      obj.deltas   = zeros( length(obj.t_total), 1 );
      obj.delta_ds = zeros( length(obj.t_total), 1 );
      obj.delta_d_ini = zeros( length(obj.t_total), 1 );

      % Build orbit in I, get contact forces, retrieve the contact angle and get
      % relative radial velocity
      for i = 1:length(obj.t_total)
        y_i = obj.y_total(i,:)';

        T_gamma = obj.s.T_gam(obj.y_total(i,1));
        T_beta  = obj.s.T_bet(obj.y_total(i,3));
        obj.pos1(:,i) = T_gamma' * (T_beta'*[0; 0; obj.s.l_OC]);

        if obj.s_total(i) == 0, state = 0; else, state = 1;  end

        [obj.F_cxs(i), obj.F_cys(i), ~, obj.delta_d_ini(i)] = ...
          contact_force(y_i, obj.s, obj.cmod, state);

        obj.alphas(i) = obj.s.contact_ang(y_i);
        obj.deltas(i) = obj.s.calc_gap(y_i);
        [ obj.delta_ds(i), obj.v_rel_r(:,i) ] = obj.s.pen_rate(y_i);
      end

      % Scale down velocity vector
      obj.v_rel_r = obj.v_rel_r*1/15;

      % Calculate resultant radial contact force
      obj.F_r = sqrt(obj.F_cxs.^2 + obj.F_cys.^2);
    end


    % 'setupPlots' sets up the different plot windows.
    function setupPlots(obj)

      % Create main figure window
      obj.mainfig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

      % Rotor and stator perimeter circles
      subplot(2,2,[1,3]); hold on

      obj.rotor  = plot(0,0,'r','LineWidth',1.5); grid on
      obj.stator = plot(0,0,'k','LineWidth',1.2);

      % Center cross markers
      obj.rotor_cross = plot(0,0,'r','LineWidth',1.5);

      obj.stator_cross = plot(obj.y_total(1,7), obj.y_total(1,9),'k+',...
        'LineWidth',1.5,'MarkerSize',20);

      % Revolution mark
      obj.rev_mark = plot(0,0,'r','LineWidth',2);

      % Quiver plot
      obj.v_vector = quiver( 0, 0, 0, 0, ...
        'LineWidth', 1.2, 'AutoScale', 'on', ...
        'MaxHeadSize', 0.9 );

      % Define axes
      axis(1*[-0.016 0.016 -0.016 0.016]);
      axis equal
      xticklabels([]); yticklabels([]);

      % Force plot
      obj.fplt = subplot(2,2,2);
      obj.force_plt = plot(0,0,'b','LineWidth',1.5); grid on

      xlabel('Time [s]')
      ylabel('Force [N]')
      obj.fplt.YLim = [0 max(obj.F_r)];

      % Force-indentation plot
      subplot(2,2,4);
      obj.forcein_plt = plot(0,0,'k','LineWidth',1.5); grid on
      xlabel('Indentation [mm]')
      ylabel('Normal force [N]')
    end


    % 'setupGUI' sets up GUI elements.
    function setupGUI(obj)
      % Add info text field
      obj.info_field = uicontrol('style', 'text', ...
        'backgroundcolor', [1 1 1], ...
        'fontname', 'Monospaced', ...
        'fontsize', 14, ...
        'horizontalalignment', 'left', ...
        'units', 'normalized', ...
        'position', [0.01, 0.58, 0.12, 0.36]);

      % Add impact counter in force-indentation plot
      obj.imp_field = uicontrol('style', 'text', ...
        'string', ['0/', num2str(ceil(length(obj.te_total)/2))], ...
        'backgroundcolor', [1 1 1], ...
        'fontname', 'Monospaced', ...
        'fontsize', 14, ...
        'horizontalalignment', 'left', ...
        'units', 'normalized', ...
        'position', [0.58, 0.42, 0.05, 0.03]);

      % Add slider control
      uicontrol('style','text', 'string', 'Time step', 'fontsize', 14, ...
        'fontname', 'Monospaced', ...
        'horizontalalignment', 'left', ...
        'units','normalized', 'position',[0.015 0.04 0.1 0.05])

      obj.sld_handle = uicontrol('style','slider',...
        'backgroundcolor', [1 1 1], ...
        'units','normalized',...
        'position',[0.01 0.04 0.4 0.01],...
        'min',1,'max',length(obj.t_total),'val',1,...
        'SliderStep', [30/length(obj.t_total) 30/length(obj.t_total)]);

      % Create eventlistener for a change of slider value
      addlistener(obj.sld_handle,'Value','PostSet', @(~,~) obj.updatefig());

      % Add slider sensitivity control
      obj.sens_button = uicontrol('style','togglebutton', ...
        'string','Fine Step','fontsize', 12, ...
        'units','normalized', ...
        'position',[0.44 0.025 0.06 0.0375]);

      % Eventlistener for sens control
      addlistener(obj.sens_button, 'Value', 'PostSet', @(~,~)obj.toggleSldSens);

      % Add buttons the jump forward and backward in impacts
      obj.fwd_button = uicontrol('style', 'pushbutton', 'string', 'Next', ...
        'fontsize', 12, ...
        'units', 'normalized', ...
        'position',[0.06 0.10 0.04 0.0375]);

      obj.bwd_button = uicontrol('style', 'pushbutton', 'string', 'Previous', ...
        'fontsize', 12, ...
        'units', 'normalized', ...
        'position',[0.015 0.10 0.04 0.0375]);

      % Eventlisteners
      set(obj.fwd_button, 'Callback', ...
        @(~, ~) obj.stepImpact(obj.fwd_button.String));
      set(obj.bwd_button, 'Callback', ...
        @(~, ~) obj.stepImpact(obj.bwd_button.String));

      % Add listener for keyboard shortcuts
      set(obj.mainfig, 'KeyPressFcn', ...
        @(src, event) obj.keypress( src, event ));
      set(obj.sld_handle, 'KeyPressFcn', ...
        @(src, event) obj.keypress( src, event ));

      % Focus the slider
      uicontrol(obj.sld_handle)

    end

    function updatefig(obj)
      % 'updatefig' updates the debug plot.

      % Get plot index from slider
      i = round(get(obj.sld_handle, 'Value'));

      % Determine new positions
      rotor_x  = obj.s.r_r*cos(linspace(0,2*pi)) + obj.pos1(1,i);
      rotor_y  = obj.s.r_r*sin(linspace(0,2*pi)) + obj.pos1(2,i);
      stator_x = obj.s.r_s*cos(linspace(0,2*pi)) + obj.y_total(i,7);
      stator_y = obj.s.r_s*sin(linspace(0,2*pi)) + obj.y_total(i,9);

      % Update rotor cross
      rotorcr_x = [obj.pos1(1,i), 7e-4*cos(obj.y_total(i,5))+obj.pos1(1,i), ...
        obj.pos1(1,i), 7e-4*cos(obj.y_total(i,5)+0.5*pi)+obj.pos1(1,i), ...
        obj.pos1(1,i), 7e-4*cos(obj.y_total(i,5)+pi    )+obj.pos1(1,i), ...
        obj.pos1(1,i), 7e-4*cos(obj.y_total(i,5)+1.5*pi)+obj.pos1(1,i)];

      rotorcr_y = [obj.pos1(2,i), 7e-4*sin(obj.y_total(i,5))+obj.pos1(2,i), ...
        obj.pos1(2,i), 7e-4*sin(obj.y_total(i,5)+0.5*pi)+obj.pos1(2,i), ...
        obj.pos1(2,i), 7e-4*sin(obj.y_total(i,5)+pi    )+obj.pos1(2,i), ...
        obj.pos1(2,i), 7e-4*sin(obj.y_total(i,5)+1.5*pi)+obj.pos1(2,i)];

      % Update revolution mark
      rev_x = [(obj.s.r_r/1.3)*cos(obj.y_total(i,5))+obj.pos1(1,i), ...
        obj.s.r_r*cos(obj.y_total(i,5))+ obj.pos1(1,i)];
      rev_y = [(obj.s.r_r/1.3)*sin(obj.y_total(i,5))+ obj.pos1(2,i), ...
        obj.s.r_r*sin(obj.y_total(i,5))+ obj.pos1(2,i)];

      % Force plot data
      if i > 200 && i < (length(obj.t_total) - 200)
        set(obj.force_plt, 'XData', obj.t_total(i-200:i), ...
                           'YData', obj.F_r(i-200:i))
        set(obj.fplt, 'XLim', [obj.t_total(i)-1e-3 obj.t_total(i)+1e-3])
      end

      % Update
      set(obj.rotor, 'XData', rotor_x,  'YData', rotor_y);
      set(obj.stator, 'XData', stator_x, 'YData', stator_y);
      set(obj.rotor_cross, 'XData', rotorcr_x, 'YData', rotorcr_y);
      set(obj.stator_cross, 'XData', obj.y_total(i,7), ...
                            'YData', obj.y_total(i,9));
      set(obj.rev_mark, 'XData', rev_x, 'YData', rev_y);

      % Debug values
      if obj.s_total(i) == 0, state = 'False'; else, state = 'True';  end

      info_str = sprintf([' DEBUG VALUES:\n\n time = %.5f s\n', ...
        ' F_cx = %7.1f N\n F_cy = %7.1f N\n F_r  = %7.1f N\n\n', ...
        ' delta   = %10.7f\n delta_d = %10.7f\n', ...
        ' deltain = %10.7f\n\n Contact: %s'], ...
        obj.t_total(i), obj.F_cxs(i), obj.F_cys(i), obj.F_r(i), ...
        obj.deltas(i), obj.delta_ds(i), obj.delta_d_ini(i), state);

      set(obj.info_field, 'string', info_str);

      % Velocity vector
      % (XSTART, YSTART) (XDIR, YDIR)
      vector = [obj.pos1(1,i) obj.pos1(2,i) obj.v_rel_r(1,i) obj.v_rel_r(2,i)];

      set(obj.v_vector,'XData',vector( :, 1 ), 'YData', vector( :, 2), ...
                       'UData',vector( :, 3 ), 'VData', vector( :, 4))

      % Force-indentation plots
      % Find number of impact events that have occured so far
      im_idx = length(find(obj.te_total <= obj.t_total(i)));

      % If the number above is uneven we're currently in contact
      if mod(im_idx, 2) == 1

        % Start and end indices for the contact state change
        start_idx = find(obj.t_total == obj.te_total(im_idx));

        % Check weather simulation ended during impact
        if length(obj.te_total) > im_idx
          end_idx = find(obj.t_total == obj.te_total(im_idx+1));
        else
          end_idx = length(obj.t_total);
        end

        % Update identation plot
        forcein_x = obj.deltas(start_idx:end_idx);
        forcein_y = obj.F_r(start_idx:end_idx);
        set(obj.forcein_plt, 'XData', forcein_x, 'YData', forcein_y)
        set(obj.imp_field, 'string', ...
          [num2str(floor(im_idx/2)+1), '/', ...
           num2str(ceil(length(obj.te_total)/2))])
      end

      drawnow
    end

    function keypress(obj, ~, event)
      % 'keypress' is a function defining shortcut actions for the debug window.

      % Make the 'f' key toggle the sensitivity button
      if strcmp(event.Key, 'f')
        if obj.sens_button.Value == 1
          set(obj.sens_button, 'Value', 0)
        else
          set(obj.sens_button, 'Value', 1)
        end
      elseif strcmp(event.Key, 's')
        % Reset slider to start
        obj.sld_handle.Value = 1;
      elseif strcmp(event.Key, 'n')
        obj.stepImpact('Next')
      elseif strcmp(event.Key, 'p')
        obj.stepImpact('Previous')
      end
    end

    function toggleSldSens(obj)
      % 'toggleSldSens' is a callback function for a sens_button.  It toggles
      % the slider step.

      if obj.sens_button.Value == 1
        set(obj.sld_handle, 'SliderStep', ...
          [1/length(obj.t_total) 30/length(obj.t_total)]);
      else
        set(obj.sld_handle, 'SliderStep', ...
          [30/length(obj.t_total) 30/length(obj.t_total)]);
      end

      % Focus slider
      uicontrol(obj.sld_handle)
    end

    function stepImpact(obj, dir)
      % 'stepImpact' is a callback function to jump to the next/previous
      % impacts in a debug gui.

      % Get plot index from slider
      i = round(get(obj.sld_handle, 'Value'));

      % Find number of impact events that have occured so far
      im_idx = length(find(obj.te_total <= obj.t_total(i)));

      % Determine steps, in terms of number of events, to increase/decrease
      if mod(im_idx, 2) == 0 || im_idx == 0
        step_idx = 1;
      else
        step_idx = 2;
      end

      % Jump forward or backward
      if strcmp(dir, 'Next')
        % If not at the end
        if length(obj.te_total) > im_idx + step_idx
          new_value = find( obj.t_total == obj.te_total(im_idx + step_idx) );
        else
          uicontrol(obj.sld_handle)
          return
        end
      elseif strcmp(dir, 'Previous')
        % If not at the beginning
        if im_idx >= 2
          new_value = find( obj.t_total == obj.te_total(im_idx - step_idx) );
        else
          uicontrol(obj.sld_handle)
          return
        end
      end

      % Update slider
      obj.sld_handle.Value = new_value;

      % Refocus the slider
      uicontrol(obj.sld_handle)
    end
  end % methods
end % class

