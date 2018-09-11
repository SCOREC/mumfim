function convert_network(old_filename, new_filename)
%% extract from old file
old_file = dlmread(old_filename,' ');
 
nnd = old_file(1,1);
old_file(1,:) = [];
old_file(:,1) = [];
old_file(all(old_file==0,2),:) = [];
 
nnds = zeros(nnd,3);
new_edges = old_file(:,1:2) - 1;
 
nds1 = old_file(:,[1 3 4 5]);
nds1 = unique(nds1,'rows');
nds2 = old_file(:,[2 6 7 8]);
nds2 = unique(nds2,'rows');
nnds(nds1(:,1),:) = nds1(:,2:4);
nnds(nds2(:,1),:) = nds2(:,2:4);
 
 
%% print out to new file
new_file = fopen(new_filename,'w');
nodes = size(nnds,1);
edges = size(new_edges,1);
periodics = 0;
fprintf(new_file,'%i %i %i\n', nodes, edges, periodics);
for n = 1 : nodes
   fprintf(new_file,'%f %f %f\n', nnds(n,1), nnds(n,2), nnds(n,3));
end
for n = 1 : edges
   fprintf(new_file,'%i %i\n', new_edges(n,1), new_edges(n,2));
end
%fprintf(new_file,'%f\n',0.5);
fclose(new_file);
end
