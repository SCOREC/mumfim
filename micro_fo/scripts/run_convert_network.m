function run_convert_network(directory)
  % converts networks and generate parameter files
  % INPUTS:
  %   directory - directory where the networks to convert reside

  files = dir(fullfile(directory, '*.txt'));

  names = {files.name};
  for i=1:size(names,2)
    names{i}
    % append new to the name before the numeral
    underscores = strfind(names{i},'_');
    lstUnd = underscores(end);
    new_name = strcat(directory, '/' ,names{i}(1:lstUnd), 'new', names{i}(lstUnd:end))
    convert_network(names{i}, new_name);
    generate_parameters(new_name, 5E-5, 1E8);
  end
end
