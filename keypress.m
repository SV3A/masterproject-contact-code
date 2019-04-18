function keypress(~, event, sens_but_handle)
  % 'keypress' is a function defining shortcut actions

  % Make the 'f' key toggle the sensitivity button
  if strcmp(event.Key, 'f')
    if sens_but_handle.Value == 1
      set(sens_but_handle, 'Value', 0)
    else
      set(sens_but_handle, 'Value', 1)
    end
  end
end

