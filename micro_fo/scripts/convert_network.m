function convert_network(old_filename, new_filename)
%% extract from old file
old_file = dlmread(old_filename,' ',1,1);

old_file(all(old_file==0,2),:)=[];
[u1,uidx,~] = unique(old_file(:,1));
new_nodes = old_file(uidx,3:5);

[u2,uidx2,~] = unique(old_file(:,2));
uidx2(ismember(u2,u1))=[];
new_nodes = [new_nodes ; old_file(uidx2,6:8)];

new_edges = old_file(:,1:2) - 1;

%% print out to new file
new_file = fopen(new_filename,'w');
nodes = size(new_nodes,1);
edges = size(new_edges,1);
periodics = 0;
fprintf(new_file,'%i %i %i\n', nodes, edges, periodics);
for n = 1 : nodes
   fprintf(new_file,'%f %f %f\n', new_nodes(n,1), new_nodes(n,2), new_nodes(n,3));
end
for n = 1 : edges
   fprintf(new_file,'%i %i\n', new_edges(n,1), new_edges(n,2));
end
fprintf(new_file,'%f\n',0.5);
fclose(new_file);
end