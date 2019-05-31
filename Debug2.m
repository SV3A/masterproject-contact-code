classdef Debug2 < handle
  % 'Debug1' is a class used to create a visual debugging tool.

  properties (Access = private)
    % Solution data

    time           % time
    rot_x          % x-coord of rotor in contact plane
    rot_y          % y-coord of rotor in contact plane
    sta_x          % x-coord of stator in contact plane
    sta_y          % y-coord of stator in contact plane
    theta          % Rotation angle
    fn             % Magnitude of normal force
    d              % Indentation
    r_r = 25e-3/2; % Rotor radius [m]
    r_s = 29e-3/2; % Stator radius [m]


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
    sld_step
    sld_finestep
  end % properties


  methods (Access = public)

    function obj = Debug2(time, rot_x, rot_y, sta_x, sta_y, theta, fn, d)
      % Constructor function.

      obj.time  = time;
      obj.rot_x = rot_x;
      obj.rot_y = rot_y;
      obj.sta_x = sta_x;
      obj.sta_y = sta_y;
      obj.theta = theta;
      obj.fn    = fn;
      obj.d     = d;

      % Set slider step
      obj.sld_step = 30/length(obj.time);
      obj.sld_finestep = 1/length(obj.time);

      % Setup debug plot and start GUI
      obj.setupPlots();
      obj.setupGUI();
    end

  end % methods


  methods (Access = private)

    function setupPlots(obj)
      % 'setupPlots' sets up the different plot windows.

      % Create main figure window
      obj.mainfig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

      % Rotor and stator perimeter circles
      subplot(2,2,[1,3]); hold on

      obj.rotor  = plot(0,0,'r','LineWidth',1.5); grid on
      obj.stator = plot(0,0,'k','LineWidth',1.2);

      % Center cross markers
      obj.rotor_cross = plot(0,0,'r','LineWidth',1.5);

      obj.stator_cross = plot(obj.rot_x(1), obj.rot_y(1),'k+',...
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
      obj.fplt.YLim = [0 max(obj.fn)];

      % Force-indentation plot
      subplot(2,2,4);
      obj.forcein_plt = plot(0,0,'k','LineWidth',1.5); grid on
      xlabel('Indentation [mm]')
      ylabel('Normal force [N]')
    end


    function setupGUI(obj)
      % 'setupGUI' sets up GUI elements.

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
        'string', ['0/'],... %num2str(ceil(length(obj.te_total)/2))], ...
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
        'min',1,'max',length(obj.time),'val',1,...
        'SliderStep', [obj.sld_step obj.sld_step]);

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

      obj.bwd_button = uicontrol('style', 'pushbutton', ...
        'string', 'Previous', 'fontsize', 12, ...
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
      rotor_x  = obj.r_r*cos(linspace(0,2*pi)) + obj.rot_x(i);
      rotor_y  = obj.r_r*sin(linspace(0,2*pi)) + obj.rot_y(i);
      stator_x = obj.r_s*cos(linspace(0,2*pi)) + obj.sta_x(i);
      stator_y = obj.r_s*sin(linspace(0,2*pi)) + obj.sta_y(i);

      % Update rotor cross
      rotorcr_x = [obj.rot_x(i), 7e-4*cos(obj.theta(i))+obj.rot_x(i), ...
        obj.rot_x(i), 7e-4*cos(obj.theta(i)+0.5*pi)+obj.rot_x(i), ...
        obj.rot_x(i), 7e-4*cos(obj.theta(i)+pi    )+obj.rot_x(i), ...
        obj.rot_x(i), 7e-4*cos(obj.theta(i)+1.5*pi)+obj.rot_x(i)];

      rotorcr_y = [obj.rot_y(i), 7e-4*sin(obj.theta(i))+obj.rot_y(i), ...
        obj.rot_y(i), 7e-4*sin(obj.theta(i)+0.5*pi)+obj.rot_y(i), ...
        obj.rot_y(i), 7e-4*sin(obj.theta(i)+pi    )+obj.rot_y(i), ...
        obj.rot_y(i), 7e-4*sin(obj.theta(i)+1.5*pi)+obj.rot_y(i)];

      % Update revolution mark
      rev_x = [(obj.r_r/1.3)*cos(obj.theta(i))+obj.rot_x(i), ...
        obj.r_r*cos(obj.theta(i))+ obj.rot_x(i)];
      rev_y = [(obj.r_r/1.3)*sin(obj.theta(i))+ obj.rot_y(i), ...
        obj.r_r*sin(obj.theta(i))+ obj.rot_y(i)];

      % Force plot data
      if i > 200 && i < (length(obj.time) - 200)
        set(obj.force_plt, 'XData', obj.time(i-200:i), ...
                           'YData', obj.fn(i-200:i))
        set(obj.fplt, 'XLim', [obj.time(i)-1e-3 obj.time(i)+1e-3])
      end

      % Update
      set(obj.rotor, 'XData', rotor_x,  'YData', rotor_y);
      set(obj.stator, 'XData', stator_x, 'YData', stator_y);
      set(obj.rotor_cross, 'XData', rotorcr_x, 'YData', rotorcr_y);
      set(obj.stator_cross, 'XData', obj.sta_x(i), ...
                            'YData', obj.sta_y(i));
      set(obj.rev_mark, 'XData', rev_x, 'YData', rev_y);

      % Debug values
      if 0 == 0, state = 'False'; else, state = 'True';  end

      lolz = 'N/A';
      info_str = sprintf([' DEBUG VALUES:\n\n time = %.5f s\n', ...
        ' F_cx = %7.1f N\n F_cy = %7.1f N\n F_r  = %7.1f N\n\n', ...
        ' delta   = %10.7f\n delta_d = %10.7f\n', ...
        ' deltain = %10.7f\n\n Contact: %s'], ...
        lolz, lolz, lolz, lolz, ...
        lolz, lolz, lolz, state);

      set(obj.info_field, 'string', info_str);

      % Velocity vector
      % (XSTART, YSTART) (XDIR, YDIR)
      %vector = [obj.rot_x(i) obj.rot_y(i) obj.v_rel_r(1,i) obj.v_rel_r(2,i)];

      %set(obj.v_vector,'XData',vector( :, 1 ), 'YData', vector( :, 2), ...
                       %'UData',vector( :, 3 ), 'VData', vector( :, 4))

      % Force-indentation plots
      % Find number of impact events that have occured so far
      %im_idx = length(find(obj.te_total <= obj.time(i)));

      %% If the number above is uneven we're currently in contact
      %if mod(im_idx, 2) == 1

        %% Start and end indices for the contact state change
        %start_idx = find(obj.time == obj.te_total(im_idx));

        %% Check weather simulation ended during impact
        %if length(obj.te_total) > im_idx
          %end_idx = find(obj.time == obj.te_total(im_idx+1));
        %else
          %end_idx = length(obj.time);
        %end

        %% Update identation plot
        %forcein_x = obj.deltas(start_idx:end_idx);
        %forcein_y = obj.F_r(start_idx:end_idx);
        %set(obj.forcein_plt, 'XData', forcein_x, 'YData', forcein_y)
        %set(obj.imp_field, 'string', ...
          %[num2str(floor(im_idx/2)+1), '/', ...
           %num2str(ceil(length(obj.te_total)/2))])
      %end

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
        set(obj.sld_handle, 'SliderStep', [obj.sld_finestep obj.sld_step]);
      else
        set(obj.sld_handle, 'SliderStep', [obj.sld_step obj.sld_step]);
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
      im_idx = length(find(obj.te_total <= obj.time(i)));

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
          new_value = find( obj.time == obj.te_total(im_idx + step_idx) );
        else
          uicontrol(obj.sld_handle)
          return
        end
      elseif strcmp(dir, 'Previous')
        % If not at the beginning
        if im_idx >= 2
          new_value = find( obj.time == obj.te_total(im_idx - step_idx) );
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
