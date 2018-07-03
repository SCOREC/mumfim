clear, clc;
% converts networks and generate parameter files
d = '/fasttmp/deoges/develop/bio_probs/fiber_networks/set_4/';
files = dir(fullfile(d, '*.txt'));

names = {files.name};
for i=1:size(names,2)
  % append new to the name before the numeral
  underscores = strfind(names{i},'_');
  lstUnd = underscores(end);
  new_name = strcat(names{i}(1:lstUnd), 'new', names{i}(lstUnd:end))
  convert_network(names{i}, new_name);
  generate_parameters(new_name, 5E-5, 1E8);
end
