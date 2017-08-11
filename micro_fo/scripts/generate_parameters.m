function generate_parameters(filename)
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
fprintf(param_file,'%i %e %f\n', 0, 3.49911271e-08,10000); % this is the reaction type and params
for n = 1 : edges
   fprintf(param_file,'%i\n', 0); % this is the index of the reaction, not the reaction TYPE
end
fclose(param_file);
end