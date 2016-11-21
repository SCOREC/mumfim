function plot_supported_networks(fname, num_networks_start, num_networks_end, color)

for i=num_networks_start:num_networks_end
    input_file = [fname '_' num2str(i) '.txt'];
    [nodes,fibers,CLIP_BOUNDARY] = read_supported_network( input_file );
    fibers(:,2:3) = fibers(:,2:3)+1;% temporary fix, should output element number from 1 rather than 0.
    plot_supported_mesh(nodes(:,2:4),fibers(:,2:4),i+1,color);
end
end



function [ nodes, fibers, clip_boundary ] = read_supported_network( filename )

% Copied from network clip function
% Except checks for optional clip boundary value included at end of file

nums = dlmread(filename, ' ', [0 0 0 2]);

raw_data = dlmread( filename, ' ', [1 0 nums(3) 9]); % load vect file skip header

% clip_boundary = dlmread( filename, ' ', [nums(3)+1 0 nums(3)+1 0]);
clip_boundary = 0.5;
fiber_number                = raw_data( : , 1 );
fiber_node_a                = raw_data( : , 2 );
fiber_node_b                = raw_data( : , 3 );

fiber_node_a_x              = raw_data( : , 4 );
fiber_node_a_y              = raw_data( : , 5 );
fiber_node_a_z              = raw_data( : , 6 );

fiber_node_b_x              = raw_data( : , 7 );
fiber_node_b_y              = raw_data( : , 8 );
fiber_node_b_z              = raw_data( : , 9 );                      

fiber_flag                  = raw_data( : , 10);

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
                                                             
fibers                      = [fiber_number fiber_node_a fiber_node_b fiber_flag];



end


function plot_supported_mesh( nodes, fibers, fign, color)

% INPUT
% =====
% NODES     N-by-3 array of xyz coordinates for N nodes
% FIBERS    N-by-3 array of start-end nodes for N fibers with fiber flag

figure(fign);

plot3( nodes(:,1) , nodes(:,2), nodes(:,3),'o', ...  
    'LineWidth',1, 'MarkerEdgeColor',[0,0,0], ...
    'MarkerFaceColor', [0,0,0], 'MarkerSize',8 );
box on
hold on;

for n = 1 : size( fibers, 1 ) % count rows
    
    if fibers(n,3) == 0
        x(1) = nodes( fibers(n,1), 1 ); % node 1 x coord
        y(1) = nodes( fibers(n,1), 2 ); % node 1 y coord
        z(1) = nodes( fibers(n,1), 3 ); % node 1 z coord
        
        x(2) = nodes( fibers(n,2), 1 ); % node 2 x coord
        y(2) = nodes( fibers(n,2), 2 ); % node 2 y coord
        z(2) = nodes( fibers(n,2), 3 ); % node 2 z coord
        
        plot3( x, y, z, 'color',color{fibers(n,3)+1} );
    elseif fibers(n,3) == 1    
        x(1) = nodes( fibers(n,1), 1 ); % node 1 x coord
        y(1) = nodes( fibers(n,1), 2 ); % node 1 y coord
        z(1) = nodes( fibers(n,1), 3 ); % node 1 z coord
        
        x(2) = nodes( fibers(n,2), 1 ); % node 2 x coord
        y(2) = nodes( fibers(n,2), 2 ); % node 2 y coord
        z(2) = nodes( fibers(n,2), 3 ); % node 2 z coord
        
        plot3( x, y, z, 'color',color{fibers(n,3)+1} );
    end
        
end

set(gcf, 'color', 'white');

axis equal;

end

