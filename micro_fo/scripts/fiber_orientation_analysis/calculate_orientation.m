function [fmain, fsupport] = calculate_orientation(filename)

[nodes, fibers, clip_boundary] = read_network( filename );

% Initialize nodes and elements
nodes = nodes(:,2:4);

elem = fibers(:,2:4); % also read in fiber flag at end.
elem(:,1:2) = elem(:,1:2)+1;
lengths = calc_lengths(nodes,elem);
cosine_main = [];
cosine_support = [];
cosine = zeros(size(elem,1),1);
tension_axis = [1,0,0];
main_idx = 1; support_idx = 1;
for i=1:size(elem,1)

        n1 = elem(i,1);
        n2 = elem(i,2);

        cx = (nodes(n2,1) - nodes(n1,1)) / lengths(i);
        cy = (nodes(n2,2) - nodes(n1,2)) / lengths(i);
        cz = (nodes(n2,3) - nodes(n1,3)) / lengths(i);
        
        cosine(i) = cx*tension_axis(1) + cy*tension_axis(2) + cz*tension_axis(3);
        if elem(i,3) == 0
            cosine_main(main_idx) = cosine(i);
            main_idx = main_idx + 1;
        elseif elem(i,3) == 1
            cosine_support(support_idx) = cosine(i);
            support_idx = support_idx + 1;
        end
end  

% Some values in the orientations array are less than 0. 
% These values that are less than 0 represent fibers that are more than 
% 180 degrees from the tension axis. 

% When calculating the orientation factor, we are concerned with the square
% of the values of the orientations array. Therefore, the direction of the
% fibers are not taken into account.

% The orientation factor, f, is described in Eq (3.59) of Dynamics of Fibre
% Formation and Processing by Roland Beyreuther and Harald Brunig.
% f = 0: fibers are randomly oriented.
% f = 1: fibers are aligned parallel to the tension axis.
cos2main = mean(cosine_main.^2);
cos2support = mean(cosine_support.^2);
fmain = (3*cos2main - 1)/2;
fsupport = (3*cos2support - 1)/2;
end

function [ nodes, fibers, clip_boundary ] = read_network( filename )

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

function [lengths] = calc_lengths(nodes,elem)

lengths = zeros(size(elem,1),1);
for i=1:size(elem,1)
    lengths(i) = sqrt( (nodes(elem(i,2),1)-nodes(elem(i,1),1))^2 + ...
                       (nodes(elem(i,2),2)-nodes(elem(i,1),2))^2 + ... 
                       (nodes(elem(i,2),3)-nodes(elem(i,1),3))^2 );

end

end