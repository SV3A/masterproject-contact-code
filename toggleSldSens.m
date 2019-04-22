function toggleSldSens(but_value, sld_handle, t_total)
% 'toggleSldSens' is a callback function for a button.  It toggles the slider
% step.

  if but_value == 1
    set(sld_handle, 'SliderStep', [1/length(t_total) 30/length(t_total)]);
    uicontrol(sld_handle)
  else
    set(sld_handle, 'SliderStep', [30/length(t_total) 30/length(t_total)]);
    uicontrol(sld_handle)
  end
end

