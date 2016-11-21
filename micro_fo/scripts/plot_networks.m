function plot_networks(fname, num_networks)

for i=1:num_networks
    input_file = [fname '_' num2str(i) '.txt'];
    [nodes,fibers,CLIP_BOUNDARY] = read_network( input_file );
    plot_mesh(nodes(:,2:4),fibers(:,2:3),i);
end
end



function [ nodes, fibers, clip_boundary ] = read_network( filename )

% Copied from network clip function
% Except checks for optional clip boundary value included at end of file

nums = dlmread(filename, ' ', [0 0 0 2]);

raw_data = dlmread( filename, ' ', [1 0 nums(3) 8]); % load vect file skip header

clip_boundary = dlmread( filename, ' ', [nums(3)+1 0 nums(3)+1 0]);

fiber_number                = raw_data( : , 1 );
fiber_node_a                = raw_data( : , 2 );
fiber_node_b                = raw_data( : , 3 );

fiber_node_a_x              = raw_data( : , 4 );
fiber_node_a_y              = raw_data( : , 5 );
fiber_node_a_z              = raw_data( : , 6 );

fiber_node_b_x              = raw_data( : , 7 );
fiber_node_b_y              = raw_data( : , 8 );
fiber_node_b_z              = raw_data( : , 9 );                      

node_stack                  = [ fiber_node_a ...   % column 1 node number
                                fiber_node_a_x ... % column 2 node x coord
                                fiber_node_a_y ... % column 3 node y coord
                                fiber_node_a_z;    % column 4 node z coord
                                
                                fiber_node_b ...   % node b piled
                                fiber_node_b_x ... % underneath the
                                fiber_node_b_y ... % first node
                                fiber_node_b_z ];
                            
nodes                       = unique( node_stack , 'rows' ); % unique rows
                                                             % ordered from
                                                             % least to 
                                                             % greatest via
                                                             % the unique()
                                                             % function
                                                             
fibers                      = [fiber_number fiber_node_a fiber_node_b];



end


function plot_mesh( nodes, fibers, fign)

% INPUT
% =====
% NODES     N-by-3 array of xyz coordinates for N nodes
% FIBERS    N-by-2 array of start-end nodes for N fibers

figure(fign);

plot3( nodes(:,1) , nodes(:,2), nodes(:,3),'o', ...  
    'LineWidth',1, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor', [0.5 0.5 0.8], 'MarkerSize',3 );

hold on;

for n = 1 : size( fibers, 1 ) % count rows
    
    x(1) = nodes( fibers(n,1), 1 ); % node 1 x coord
    y(1) = nodes( fibers(n,1), 2 ); % node 1 y coord
    z(1) = nodes( fibers(n,1), 3 ); % node 1 z coord
        
    x(2) = nodes( fibers(n,2), 1 ); % node 2 x coord
    y(2) = nodes( fibers(n,2), 2 ); % node 2 y coord
    z(2) = nodes( fibers(n,2), 3 ); % node 2 z coord
        
    plot3( x, y, z );
        
    hold on;
        
end

set(gcf, 'color', 'white');

axis equal;

end

