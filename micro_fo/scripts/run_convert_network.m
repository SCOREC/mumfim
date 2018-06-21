clear, clc;
% converts networks and generate parameter files
d = '/fasttmp/mersoj/develop/biotissue_problems/fiber_networks/density_100_networks/';
files = dir(fullfile(d, '*.txt'));

names = {files.name}
for i=1:size(names,2)
  names{i}
  % append new to the name before the numeral
  underscores = strfind(names{1},'_');
  lstUnd = underscores(end);
  new_name = strcat(names{i}(1:lstUnd), 'new', names{i}(lstUnd:end));
  convert_network(names{i}, new_name);
  generate_parameters(new_name, 3.846494e-15, 10000);
end
