function export_values(par_list, values)
% 'export_value' writes values to a text file.
% INPUT:
%   par_list    : Parameter list e.g. par_list = {'t', 'rotor_x', ...}
%   value_vector: All values to be printed in one long 1xn vector

  % Open file
  fileID = fopen('exp.txt', 'w');

  % Format spec
  header_format = [strtrim(repmat('%-9s ' , 1, size(par_list, 2))), '\n'];
  number_format = [strtrim(repmat('%9.7f ', 1, size(par_list, 2))), '\n'];

  % Write header
  fprintf(fileID, header_format, par_list{:});

  % Write values
  fprintf(fileID, number_format, values);

  % Close file
  fclose(fileID);
end
