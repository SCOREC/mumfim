function generate_parameters(filename, Rf, Ef)
% generate_parameters: writes fiber reactions for new biotissue format
% INPUTS:
%   filename -- name of file with fiber network
%   Rf -- fiber radius (double)
%   Ef -- fiber youngs modulus (double)
% OUTPUT:
%   writes filename.params
% TODO: take a vector Rf and Ef (length of number of reactions),
%       reactions (length of number of fibers) that gives reaction number for each fiber

%% extract from file
file = fopen(filename,'r');
hdr_txt = fgetl(file);
hdr = str2num(hdr_txt);
%% print out to params file
params_filename = strcat(filename,'.params');
param_file = fopen(params_filename,'w');
edges = hdr(2);
reactions = 1;
fprintf(param_file,'%i %i\n', reactions, edges);
fprintf(param_file,'%i %d %f\n', 0, Rf, Ef); % this is the reaction type and params
for n = 1 : edges
   fprintf(param_file,'%i \n', 0); % this is the index of the reaction, not the reaction TYPE
end
fclose(param_file);
end
