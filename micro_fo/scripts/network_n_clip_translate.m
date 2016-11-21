function [nodes, fibers] = network_n_clip_translate(nodes, fibers, clip_boundary)
%********************************************************
%
% Following descriptions are outdated, see history...
%
%
% DESCRIPTION
% ===========
%
% NETWORK_N_CLIPPER() reads in a network file and then clips the
% fibers along a specified inner box, forming new nodes along the faces
% of that inner cube that are added back into a new, clipped network.
% Translates and then outputs the network.
%
% INPUT
% =====
%
% NETWORK_FILE -- read in via read_network() is a text file with
% a header row containing #nodes #dof #fibers, followed by N lines of
% data for N fibers containing #fibernum #startnodenum #endnodenum
% #startnodexcoord ... #endnodezcoord as used in the multiscale code
%
% KEY VARIABLES
% =============
%
% NODES is N x 4 -- col1=node index, col2-4=node x y z coordinates
%
% FIBERS is M x 3 -- col1=fiber index, col2=node a col3=node b
%
% The code should be able to handle input networks that are not numbered
% sequentially as long as nodes and fibers are assigned unique indices,
% although this has not been tested.
%
% PARAMETERS
% ==========
%
% Search for keyword "SPECIFY" to find key parameters
%
% CLIP_BOUNDS must be set to the box boundaries you want to clip to
% NETWORK_FILE must contain the name of the network file
% OUTPUT_FILE is the name that we'll write the clipped network to
%
% HISTORY
% =======
%
% March 2015 - Changed to be a function that takes and returns 'fibers'
%              and 'nodes' data structures and returns them after clipping to
%              boundary. Instead of a stand alone process that reads and
%              writes to files.
%              - Dan Fovargue
%
%
% Dec 2010 -- Created -- MFH
% Mar 2011 -- Revised to translate and check for zero length fibers -- MFH
%
%********************************************************

if size(nodes,2) ~= 4
    error('Nodes structure wrong size in clip function')
end
if size(fibers,2) ~= 3
    error('Fibers structure wrong size in clip function')
end

tic

fprintf('In clipping function (clipping to %f)...\n',clip_boundary);

clip_bounds =   [  -clip_boundary clip_boundary ...
                   -clip_boundary clip_boundary ...
                   -clip_boundary clip_boundary ];
               
% Check for zero-length fibers in input file
[min_length min_fib_num max_length max_fib_num] = check_small_fibers(nodes, fibers);

if min_length ~= 0
    fprintf('Good -- no zero length fibers in input file...\n');
else
    error('Bad -- zero length fibers in input file...\n');
end

% Get next biggest value for fiber number from col1 of the fibers file
fiber_start_index = max( fibers(:,1) ) + 1;

% Get next biggest value for node number from col1 of the node file
node_start_index = max( nodes(:,1) ) + 1;

new_node_stack      = [];
new_fiber_stack     = [];

 % Find interior fibers - completely inside/on clipping boundaries
[interior_fiber_numbers,clip_nodes,clip_fibers] = find_interior_fibers(nodes,fibers,clip_bounds);

new_node_stack      = [ new_node_stack; clip_nodes ];       % add to master stack
new_fiber_stack     = [ new_fiber_stack; clip_fibers ];     % add to master stack

% Find crossing fibers - one node inside and one node outside clip boundary
[crossing_fiber_numbers,clip_nodes,clip_fibers,...
    node_start_index,fiber_start_index ] = ...
    find_crossing_fibers(nodes, fibers, clip_bounds,...
    node_start_index,fiber_start_index);

new_node_stack      = [     new_node_stack;         % add to master stack
    clip_nodes      ];

new_fiber_stack     = [     new_fiber_stack;        % add to master stack
    clip_fibers     ];

% Spanning fibers - span through clipping boundary box but no nodes inside
[spanning_fiber_numbers,clip_nodes,clip_fibers,...
    node_start_index,fiber_start_index ] = ...
    find_spanning_fibers(nodes, fibers, clip_bounds,...
    node_start_index,fiber_start_index);

new_node_stack      = [ new_node_stack; % add to master stack
    clip_nodes ];
new_fiber_stack     = [ new_fiber_stack; % add to master stack
    clip_fibers ];

% check our new network to ensure:
%
% - our fiber numbers are consistent
% - uniqueness of fibers
% - no zero length fibers

% make sure each group of fibers is mutually exclusive
all_group_fib_num = [   interior_fiber_numbers ...
    crossing_fiber_numbers ...
    spanning_fiber_numbers ];

if length( all_group_fib_num ) ~= length( unique( all_group_fib_num ) )
    error('Bad -- overlap of fiber categories...');
else
    disp('Good -- fiber categories mutually exclusive...');
end

% check to make sure number of fibers is correct

length_fib_categories       = size( all_group_fib_num,  2 ); % count col
length_master_fib_stack     = size( new_fiber_stack,    1 ); % count row

if length_fib_categories ~= length_master_fib_stack
    error('Bad -- target fibers ~= fibers created...');
else
    disp('Good -- target fibers = fibers created...');
end

% clean up clip_nodes and clip_fibers
% nodes contains redundant nodes and non-conseq node numbers
% fibers contains non-conseq fiber numbers
[final_nodes,final_fibers] = net_cleanup(new_node_stack,new_fiber_stack);

% check for zero-length fibers again in the clipped network
[min_length min_fib_num max_length max_fib_num] = check_small_fibers(final_nodes, final_fibers);

if min_length ~= 0
    fprintf('Good -- no zero length fibers in output file...\n');
else
    error('Bad -- zero length fibers in output file...\n');
end

% translate network over a fixed amount SPECIFIY DISPLACEMENT (IF ANY) TO TRANSLATE THE NETWORK
[final_nodes,final_fibers] = net_translate( final_nodes,final_fibers, [0.0 0.0 0.0] ); 

% Determine total fiber length and volume of bounding cube
cube_volume = (clip_boundary*2)^3;
tfl = 0.0;
for i=1:size(final_fibers,1)
    x1 = final_nodes(final_fibers(i,2),2:4);
    x2 = final_nodes(final_fibers(i,3),2:4);
    tfl = tfl + norm(x1-x2,2);
end
fprintf('tfl = %f , cube volume = %f , ratio = %f \n',tfl,cube_volume,cube_volume/tfl);
fprintf('num_nodes = %i , num_elements = %i \n',size(final_nodes,1),size(final_fibers,1));
fprintf('avg fiber length = %f\n',tfl/size(final_fibers,1));

fibers = final_fibers;
nodes = final_nodes;

% Plot - just for debugging, by default is skipped
if false
    figure('Color',[1 1 1]);
    plot_network( nodes, fibers, [-0.5 0.5 -0.5 0.5 -0.5 0.5] );
end

fprintf('Finished clipping. Time to clip: %f\n',toc);

end

function plot_network( nodes, fibers, network_axis )

% INPUT
% =====
% NODES - m-by-3 array of xyz coordinates for m nodes indexed 1-m
% FIBERS - m-by-2 array of start-end node indices for 1-m fibers
% NETWORK_AXIS - [xmin xmax ymin ymax zmin zmax] plotting axis vector
%
% HISTORY
% =======
% Mar 2011 -- Last mod

nodes = nodes(:,2:4);
fibers = fibers(:,2:3);

total_nodes = size( nodes, 1 ); % num rows = num nodes

total_fibers = size( fibers, 1 ); % num rows = num fibers

min_x = min( nodes(:,1) );
max_x = max( nodes(:,1) );

min_y = min( nodes(:,2) );
max_y = max( nodes(:,2) );

min_z = min( nodes(:,3) );
max_z = max( nodes(:,3) );

% PLOT NODES
if 1 == 1
    
    plot3( nodes(:,1) , ...
        nodes(:,2), ...
        nodes(:,3),'o',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor', [rand rand rand],...
        'MarkerSize', 4);
    
    axis( network_axis );
    
    hold on;
    
end

% PLOT FIBERS
for n = 1 : total_fibers
    
    x(1) = nodes( fibers(n,1), 1 ); % node 1 x coord
    y(1) = nodes( fibers(n,1), 2 ); % node 1 y coord
    z(1) = nodes( fibers(n,1), 3 ); % node 1 z coord
    
    x(2) = nodes( fibers(n,2), 1 ); % node 2 x coord
    y(2) = nodes( fibers(n,2), 2 ); % node 2 y coord
    z(2) = nodes( fibers(n,2), 3 ); % node 2 z coord
    
    plot3( x, y, z, 'linewidth', 2);
    
    hold on;
    
end

axis( network_axis );

end


function [min_length min_fib_num max_length max_fib_num] = check_small_fibers(nodes, fibers)

% INPUT
% =====
%
% FIBERS col1 = fiber number
% FIBERS col2 = start node
% FIBERS col3 = end node
%
% NODES col1 = node number
% NODES col2 = x
% NODES col3 = y
% NODES col4 = z

total_nodes         = size( nodes, 1 );     % count rows = num nodes
total_fibers        = size( fibers, 1 );    % count rows = num fibers

for n = 1 : total_fibers
    
    fiber_num       = fibers( n , 1 );
    
    node_a_num      = fibers( n , 2 );
    
    node_b_num      = fibers( n , 3 );
    
    [node_a_x node_a_y node_a_z] = get_node_coords( nodes, node_a_num );
    
    [node_b_x node_b_y node_b_z] = get_node_coords( nodes, node_b_num );
    
    delta_coords = [node_a_x node_a_y node_a_z] - [node_b_x node_b_y node_b_z];
    
    fiber_length(n) = norm(  delta_coords ); % norm across columns
    
end

[min_length min_fib_num] = min( fiber_length );

[max_length max_fib_num] = max( fiber_length );

end


function [nodes,fibers] = net_translate( nodes, fibers, translate_vector )

% Description
% ===========
%
% NET_TRANSLATE shifts the xyz coords in NODES by a certain value.
%
% TRANSLATE_VECTOR should be a row vector for translating nodal coords.
%
% FIBERS is left unused at this point.
%
% History
% =======
%
% Dec 17 2010 -- Wrote for the clipping method -- MFH

total_nodes         = size( nodes, 1 );     % count rows = num nodes
total_fibers        = size( fibers, 1 );    % count rows = num fibers

for n = 1 : total_nodes
    
    nodes(n,2:4) = nodes(n,2:4) + [translate_vector];
    
end

end


function write_mesh( filename, nodes, fibers, CLIP_BOUNDARY)

% History
% =======
%
% Dec 01 10 -- Taken from the Delaunay generator -- MFH

total_nodes         = size( nodes, 1 );     % count rows = num nodes
total_fibers        = size( fibers, 1 );    % count rows = num fibers

tot_deg_freedom     = 3 * total_nodes;      % x y z values for nodes

fileid = fopen(filename,'w'); % writes over any existing file of this name

fprintf(fileid,'%i %i %i\n', total_nodes, tot_deg_freedom, total_fibers);

for n = 1 : total_fibers
    
    fprintf(fileid,'%i %i %i %f %f %f %f %f %f\n', ...  % file line
        n, ...                                          % fiber num
        fibers(n,2), ...                                % node 1 num
        fibers(n,3), ...                                % node 2 num
        nodes( fibers(n,2), 2 ), ...              % node 1 x coord
        nodes( fibers(n,2), 3 ), ...              % node 1 y coord
        nodes( fibers(n,2), 4 ), ...              % node 1 z coord
        nodes( fibers(n,3), 2 ), ...              % node 2 x coord
        nodes( fibers(n,3), 3 ), ...              % node 2 y coord
        nodes( fibers(n,3), 4 ));                 % node 2 z coord
    
end

fprintf(fileid,'%f\n',CLIP_BOUNDARY);

fclose( fileid );

end


function [final_nodes,final_fibers] = net_cleanup( new_node_stack, new_fiber_stack)

% first sort out AND renumber nodes

new_node_stack = sortrows( unique( new_node_stack, 'rows')); % UNIQUE() should also sort in
% order via the first column of
% values (ie node number) but
% we'll re-sort just in case

number_of_nodes = size( new_node_stack, 1); % number of rows = number nodes

new_node_list = 1 : number_of_nodes;

old_node_list = new_node_stack( : , 1 );

new_node_stack( : , 1 ) = 1 : number_of_nodes; % renumber nodes

% renumber fiber start/end nodes

number_of_fibers = size( new_fiber_stack, 1 ); % count rows = num fibers

for n = 1:number_of_fibers
    
    old_node_a = new_fiber_stack( n , 2 );
    old_node_b = new_fiber_stack( n , 3 );
    
    old_node_index_a = find( old_node_list == old_node_a );
    old_node_index_b = find( old_node_list == old_node_b );
    
    new_node_a = new_node_list( old_node_index_a );
    new_node_b = new_node_list( old_node_index_b );
    
    new_fiber_stack( n , 2 ) = new_node_a;
    new_fiber_stack( n , 3 ) = new_node_b;
    
end

% renumber fibers

new_fiber_stack(:,1) = 1 : number_of_fibers; % first col = fiber number

% output resultant fibers

final_nodes     =   new_node_stack;   % nodes cleaned up
final_fibers    =   new_fiber_stack;  % fibers cleaned up

end


function [ spanning_fiber_numbers, ...
    clip_nodes, clip_fibers, node_start_index, fiber_start_index ] = ...
    find_spanning_fibers( nodes, fibers, clip_bounds, ...
    node_start_index, fiber_start_index )

% Description
% ===========
%
% FIND_SPANNING_FIBERS() finds fibers that span the entire boundary
% box without any nodes inside the box. It does so by selecting fibers
% that intersect two different faces of the boundary box.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

number_fibers = size( fibers, 1 ); % num rows of fibers = num of fibers

number_boundary_faces = 6; % boundary faces -- always 6 -- box

[boundary_plane_pts] = gather_plane_points( clip_bounds ); % each 4 rows
% contain new
% pt coords for
% all 6 faces

spanning_fiber_numbers = []; % setup a fillable array for interior nums

clip_nodes      = []; % initialize a fillable array

clip_fibers     = []; % initialize a fillable array

for n = 1 : number_fibers % evaluate each and every fiber
    
    num_intersects = 0; % set counter for the times fiber hits face
    
    node_a_num = fibers( n , 2 );
    node_b_num = fibers( n , 3 );
    
    [node_a_x node_a_y node_a_z] = get_node_coords( nodes, node_a_num );
    [node_b_x node_b_y node_b_z] = get_node_coords( nodes, node_b_num );
    
    node_a_coords = [node_a_x node_a_y node_a_z]; % condense coords
    node_b_coords = [node_b_x node_b_y node_b_z]; % condense coords
    
    % node_a_inside is (1) if inside/on boundary (0) if outside boundary
    
    node_a_inside = check_inside_boundary( node_a_x, node_a_y, ...
        node_a_z, clip_bounds );
    
    node_b_inside = check_inside_boundary( node_b_x, node_b_y, ...
        node_b_z, clip_bounds);
    
    % initialize new_node_list for EACH fiber
    
    new_node_list       = [];
    new_node_num_list   = [];
    
    % check fiber against ALL SIX faces of the bounding cube
    
    if ( ~node_a_inside && ~node_b_inside ) % if nodes NOT inside/on box
        
        for m = 1 : number_boundary_faces % find any face intersects
            
            plane_pt_a = boundary_plane_pts( m*4 - 3, : ); % [x y z]
            plane_pt_b = boundary_plane_pts( m*4 - 2, : ); % [x y z]
            plane_pt_c = boundary_plane_pts( m*4 - 1, : ); % [x y z]
            plane_pt_d = boundary_plane_pts( m*4 - 0, : ); % [x y z]
            
            [outcome, xyz] = find_intersect( node_a_coords,...  % check to
                node_b_coords,...  % see if
                plane_pt_a,...     % fiber
                plane_pt_b,...     % hits
                plane_pt_c,...     % a face
                plane_pt_d,...
                clip_bounds,...
                m );
            
            if outcome == 1 % start tracking new

                % NEW_NODE_LIST will be ideally 2 rows 4 columns
                % after getting 2 face intersects
                % column 1 new node number
                % column 2 x coord
                % column 3 y coord
                % column 4 z coord
                
                new_node_list       = [ new_node_list;
                    node_start_index xyz];
                
                node_start_index = node_start_index + 1; % for next node
                
            end
            
            num_intersects = num_intersects + outcome;
            
        end
        
    end
    
    if num_intersects == 1
        error('Big big error -- only 1 intersect for a spanning fiber -- fix now!');
    end
    
    if num_intersects == 2 % fiber has hit 2 faces -- spans box
        
        % here we can deal with clipping the spanning fiber
        
        current_fib_num = fibers( n, 1 ); % col 1 = fiber num
        
        spanning_fiber_numbers = [spanning_fiber_numbers current_fib_num];
        
        % populate clipped output arrays
        
        clip_nodes          = [ clip_nodes;
            new_node_list ];
        
        new_node_a          = new_node_list(1,1); % node num should be here
        new_node_b          = new_node_list(2,1); % node num should be here
        
        clip_fibers         = [ clip_fibers;
            fiber_start_index new_node_a new_node_b];
        
        fiber_start_index = fiber_start_index + 1;
        
    end

end

end


function [outcome, intersect_xyz] = find_intersect(   node_a_coords,...
    node_b_coords,...
    plane_pt_a,...
    plane_pt_b,...
    plane_pt_c,...
    plane_pt_d,...
    clip_bounds,...
    face_number)

% Description
% ===========
%
% FIND_INTERSECT() locates the point at which a line segment intersects
% the face of a bounding cube. Returns OUTCOME which indicates whether
% the line intersects face and XYZ which provides the intersect point.
%
% See Wikipedia article on paramateric line/plane intersection--
% followed the parametric strategy.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

xmin = clip_bounds(1); % unpack clip_bounds array
xmax = clip_bounds(2);
ymin = clip_bounds(3);
ymax = clip_bounds(4);
zmin = clip_bounds(5);
zmax = clip_bounds(6);

x_a = node_a_coords(1); % xyz for node a on the line segment
y_a = node_a_coords(2);
z_a = node_a_coords(3);

x_b = node_b_coords(1); % xyz for node b on the line segment
y_b = node_b_coords(2);
z_b = node_b_coords(3);

x_0 = plane_pt_a(1); % xyz for one point on the plane
y_0 = plane_pt_a(2);
z_0 = plane_pt_a(3);

x_1 = plane_pt_b(1); % xyz for second point on the plane
y_1 = plane_pt_b(2);
z_1 = plane_pt_b(3);

x_2 = plane_pt_c(1); % xyz for the third point on the plane
y_2 = plane_pt_c(2);
z_2 = plane_pt_c(3);

% Solve Ax = B linear problem for the intersection of a line
% with a plane

A_matrix = [ x_a - x_b, x_1 - x_0, x_2 - x_0;
    y_a - y_b, y_1 - y_0, y_2 - y_0;
    z_a - z_b, z_1 - z_0, z_2 - z_0 ];

if det( A_matrix ) ~= 0 % do this only if A is not singular
    
    B_matrix = [   x_a - x_0;
        y_a - y_0;
        z_a - z_0 ];
    
    try
        x = A_matrix \ B_matrix; % x is our solution vector
    catch
        error('Big big error -- problem solving Ax = B in FIND_INTERSECT() -- fix now!');
    end
    
    param_t = x(1); % for parameterized line exn
    
    if (param_t >= 0 && param_t <= 1) % intersect falls btwn nodes
        % otherwise you have a line
        % outside box pointed at face
        
        intersect_xyz = node_a_coords + (node_b_coords - ...
            node_a_coords) * param_t; % find xyz of intersection
        
        intersect_x = intersect_xyz(1);
        intersect_y = intersect_xyz(2);
        intersect_z = intersect_xyz(3);
        
        % fuzz_factor to make sure xyz is near face
        
        fuzz_factor = 0.00001; % <--------------------- SELECT A FUZZ FACTOR TO SEE IF INTERSECTION POINT IS ON FACE -- MAY WANT TO DECREASE THIS VALUE
        
        % clean up intersect point -- snap firmly to a face
        % point may be a bit away from the face due to numerical
        % error from solving Ax = B above
        
        if     face_number == 1 && abs(intersect_x-xmin) < fuzz_factor
            
            intersect_x = xmin;
            
        elseif face_number == 2 && abs(intersect_x-xmax) < fuzz_factor
            
            intersect_x = xmax;
            
        elseif face_number == 3 && abs(intersect_y-ymin) < fuzz_factor
            
            intersect_y = ymin;
            
        elseif face_number == 4 && abs(intersect_y-ymax) < fuzz_factor
            
            intersect_y = ymax;
            
        elseif face_number == 5 && abs(intersect_z-zmin) < fuzz_factor
            
            intersect_z = zmin;
            
        elseif face_number == 6 && abs(intersect_z-zmax) < fuzz_factor
            
            intersect_z = zmax;
            
        end
        
        % check if intersect is somewhere ON the box
        
        if check_on_boundary( intersect_x, intersect_y, ...
                intersect_z, clip_bounds )
            
            outcome = 1; % on cube - good outcome
            
            intersect_xyz = intersect_xyz; % just a reminder of what we're trying to do here...
            
        else
            
            outcome = 0; % not on a cube face - bad outcome
            
            intersect_xyz = []; % empty out or incorrect solution
            
        end
        
    else
        
        outcome = 0; % intersect not on line btwn node a & b - bad
        
        intersect_xyz = []; % return empty set
        
    end
    
else
    
    outcome = 0; % A is a singular matrix - bad outcome
    
    intersect_xyz = []; % singular matrix - no unique intersect
    
end

end


function [boundary_plane_pts] = gather_plane_points( clip_bounds )

% Description
% ===========
%
% GATHER_PLANE_POINTS() puts together a set of points that define
% each of the faces of the boundary box that you're clipping to.
% Four points per face are used.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

xmin = clip_bounds(1); % pull out values from clip_bounds
xmax = clip_bounds(2);
ymin = clip_bounds(3);
ymax = clip_bounds(4);
zmin = clip_bounds(5);
zmax = clip_bounds(6);

xmin_face_pt_a = [xmin ymin zmin]; % put together four points on each plane
xmin_face_pt_b = [xmin ymin zmax]; % which is a face of our bounding
xmin_face_pt_c = [xmin ymax zmin]; % cube
xmin_face_pt_d = [xmin ymax zmax];

xmax_face_pt_a = [xmax ymin zmin];
xmax_face_pt_b = [xmax ymin zmax];
xmax_face_pt_c = [xmax ymax zmin];
xmax_face_pt_d = [xmax ymax zmax];

ymin_face_pt_a = [xmin ymin zmin];
ymin_face_pt_b = [xmax ymin zmin];
ymin_face_pt_c = [xmin ymin zmax];
ymin_face_pt_d = [xmax ymin zmax];

ymax_face_pt_a = [xmin ymax zmin];
ymax_face_pt_b = [xmax ymax zmin];
ymax_face_pt_c = [xmin ymax zmax];
ymax_face_pt_d = [xmax ymax zmax];

zmin_face_pt_a = [xmin ymin zmin];
zmin_face_pt_b = [xmin ymax zmin];
zmin_face_pt_c = [xmax ymax zmin];
zmin_face_pt_d = [xmax ymin zmin];

zmax_face_pt_a = [xmin ymin zmax];
zmax_face_pt_b = [xmin ymax zmax];
zmax_face_pt_c = [xmax ymax zmax];
zmax_face_pt_d = [xmax ymin zmax];

boundary_plane_pts = [ xmin_face_pt_a;
    xmin_face_pt_b;
    xmin_face_pt_c;
    xmin_face_pt_d;
    
    xmax_face_pt_a;
    xmax_face_pt_b;
    xmax_face_pt_c;
    xmax_face_pt_d;
    
    ymin_face_pt_a;
    ymin_face_pt_b;
    ymin_face_pt_c;
    ymin_face_pt_d;
    
    ymax_face_pt_a;
    ymax_face_pt_b;
    ymax_face_pt_c;
    ymax_face_pt_d;
    
    zmin_face_pt_a;
    zmin_face_pt_b;
    zmin_face_pt_c;
    zmin_face_pt_d;
    
    zmax_face_pt_a;
    zmax_face_pt_b;
    zmax_face_pt_c;
    zmax_face_pt_d  ];

end


function [ crossing_fiber_numbers, ...
    clip_nodes, clip_fibers, node_start_index, fiber_start_index ] = ...
    find_crossing_fibers( nodes, ...
    fibers, clip_bounds, ...
    node_start_index, fiber_start_index )

% Description
% ===========
%
% FIND_CROSSING_FIBERS() finds fibers that cross through ONLY one of the
% boundaries of the bounding box. Ignores fibers that have one node ON
% boundary and one node outside the boundary.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH
% Dec 01 2010 -- Changed to remove fibers with 1 node on the boundary
%                and 1 node outside the bounding box -- MFH
%

number_fibers = size( fibers, 1 ); % num rows of fibers = num of fibers

crossing_fiber_numbers = []; % setup a fillable array for interior nums

number_boundary_faces = 6; % boundary faces -- always 6 -- box

clip_nodes      = []; % initialize as an empty set
clip_fibers     = []; % initialize as an empty set

[boundary_plane_pts] = gather_plane_points( clip_bounds ); % each 4 rows
% contain new
% pt coords for
% all 6 faces

for n = 1 : number_fibers % step through each row of fibers
    
    node_a_num = fibers( n , 2 );
    node_b_num = fibers( n , 3 );
    
    [node_a_x node_a_y node_a_z] = get_node_coords( nodes, node_a_num );
    [node_b_x node_b_y node_b_z] = get_node_coords( nodes, node_b_num );
    
    node_a_coords = [node_a_x node_a_y node_a_z]; % condense coords
    node_b_coords = [node_b_x node_b_y node_b_z]; % condense coords
    
    current_fib_num = fibers(n,1);
    
    % node_a_inside is (1) if inside boundary (0) if outside or on boundary
    
    node_a_inside = check_inside_not_on_boundary( node_a_x, node_a_y, ...
        node_a_z, clip_bounds );
    
    node_b_inside = check_inside_not_on_boundary( node_b_x, node_b_y, ...
        node_b_z, clip_bounds);
    
    % node_a_bound is (1) if on boundary (0) if not on boundary
    
    node_a_bound = check_on_boundary( node_a_x, node_a_y, ...
        node_a_z, clip_bounds );
    
    node_b_bound = check_on_boundary( node_b_x, node_b_y, ...
        node_b_z, clip_bounds );
    
    if ( ~node_a_bound && ~node_b_bound ) % if nodes NOT on boundary
        
        if xor(node_a_inside,node_b_inside)  % if fiber cuts ONE boundary
            
            % we have made sure we have the fiber we want
            % now we clip the fiber--creating a new node on the boundary

            % col 1 is the fiber num - grow the list of interior fibers
            crossing_fiber_numbers = [crossing_fiber_numbers ...
                current_fib_num];
            
            % assign inside/outside roles
            if node_a_inside
                
                inside_node         = node_a_num;
                outside_node        = node_b_num;
                
                inside_node_xyz     = node_a_coords;
                outside_node_xyz    = node_b_coords;
                
            elseif node_b_inside
                
                inside_node         = node_b_num;
                outside_node        = node_a_num;
                
                inside_node_xyz     = node_b_coords;
                outside_node_xyz    = node_a_coords;
                
            end
            
            for m = 1 : number_boundary_faces % faces of clip box
                
                plane_pt_a = boundary_plane_pts(m*4-3,:); % face m pt 1 xyz
                plane_pt_b = boundary_plane_pts(m*4-2,:); % face m pt 2 xyz
                plane_pt_c = boundary_plane_pts(m*4-1,:); % face m pt 3 xyz
                plane_pt_d = boundary_plane_pts(m*4-0,:); % face m pt 4 xyz
                
                [outcome, xyz] = find_intersect( node_a_coords,...
                    node_b_coords,...  % check to see if
                    plane_pt_a,...     % fiber
                    plane_pt_b,...     % hits a
                    plane_pt_c,...     % face
                    plane_pt_d,...
                    clip_bounds,...
                    m );
                
                if outcome == 1
                    
                    break % stop trying to find intersection pts
                    
                end
                
            end
            
            new_fiber_number = fiber_start_index; % new fib num
            fiber_start_index = fiber_start_index + 1; % for next new fib
            
            new_node_number = node_start_index; % new node num
            node_start_index = node_start_index + 1; % for next new node
            
            old_node_number = inside_node;
            
            new_node_xyz    = xyz;
            old_node_xyz    = inside_node_xyz;
            
            clip_nodes      = [ clip_nodes;
                old_node_number old_node_xyz;
                new_node_number new_node_xyz ];
            
            clip_fibers     = [ clip_fibers;
                new_fiber_number ...
                old_node_number new_node_number ];
            
        end
        
    end
    
end

end


function [ interior_fiber_numbers, clip_nodes, clip_fibers ] = ...
    find_interior_fibers( nodes, fibers, clip_bounds )

% Description
% ===========
%
% FIND_INTERIOR_FIBERS() finds fibers that are contained within or on the
% clip boundaries. Returns the fiber numbers for interior fibers, not
% the row index, but the value from col1 of fibers, the unique fiber label.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

number_fibers = size( fibers, 1 ); % num rows of fibers = num of fibers

interior_fiber_numbers = []; % setup a fillable array for interior nums

clip_nodes = []; % setup a fillable array for new nodes

clip_fibers = []; % setup a fillable array for new fibers

for n = 1 : number_fibers % step through each row of fibers
    
    node_a_num = fibers( n , 2 ); % cols 2 and 3 have start/end node nums
    node_b_num = fibers( n , 3 ); % within the FIBERS array
    
    current_fib_num = fibers(n,1); % col 1 of FIBERS is the fiber num
    
    [node_a_x node_a_y node_a_z] = get_node_coords( nodes, node_a_num );
    [node_b_x node_b_y node_b_z] = get_node_coords( nodes, node_b_num );
    
    % node_a_inout is (1) if inside/on boundary (0) if outside
    
    node_a_inside = check_inside_boundary( node_a_x, node_a_y, node_a_z,...
        clip_bounds );
    
    node_b_inside = check_inside_boundary( node_b_x, node_b_y, node_b_z,...
        clip_bounds);

    if (node_a_inside && node_b_inside)  % check to see if fiber is inside
        
        % grow the list of interior fibers
        interior_fiber_numbers = [interior_fiber_numbers current_fib_num];
        
        % add to nodes and fibers
        clip_nodes =  [  clip_nodes;
            node_a_num node_a_x node_a_y node_a_z;
            node_b_num node_b_x node_b_y node_b_z ];
        
        clip_fibers = [  clip_fibers;
            current_fib_num node_a_num node_b_num ];
        
    end
    
end

end


function [node_inside] = check_inside_boundary( x, y, z, clip_bounds )

% Description
% ===========
%
% CHECK_INSIDE_BOUNDARY() checks whether a node is inside or outside the
% clipping boundary box. Returns (1) if it is (0) if it is not.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

xmin = clip_bounds(1);
xmax = clip_bounds(2);
ymin = clip_bounds(3);
ymax = clip_bounds(4);
zmin = clip_bounds(5);
zmax = clip_bounds(6);

if      x >= xmin && x <= xmax && ... % check to see if fiber
        y >= ymin && y <= ymax && ... % is contained in bounds
        z >= zmin && z <= zmax
    
    node_inside = 1; % node inside boundary
    
else
    
    node_inside = 0; % node outside OR on boundary
    
end

end


function [node_inside] = check_inside_not_on_boundary( x, y, z, ...
    clip_bounds )

% Description
% ===========
%
% CHECK_INSIDE_NOT_ON_BOUNDARY() checks whether a node is completely inside
% and not on the clipping boundary box. Returns (1) if it is (0) if it
% is not.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

xmin = clip_bounds(1);
xmax = clip_bounds(2);
ymin = clip_bounds(3);
ymax = clip_bounds(4);
zmin = clip_bounds(5);
zmax = clip_bounds(6);

if      ( (x>xmin) && (x<xmax) && ... % check to see if fiber
        (y>ymin) && (y<ymax) && ... % is contained in bounds
        (z>zmin) && (z<zmax) )
    
    node_inside = 1; % node inside boundary
    
else
    
    node_inside = 0; % node outside OR on boundary
    
end

end


function [check_result] = check_on_boundary( x, y, z, clip_bounds )

% Description
% ===========
%
% CHECK_ON_BOUNDARY() checks whether a node is on the
% clipping boundary box.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

inside_or_on        = check_inside_boundary(x,y,z,clip_bounds);
inside_not_on       = check_inside_not_on_boundary(x,y,z,clip_bounds);

if ( inside_or_on && ~inside_not_on )    % if node is inside or on the
    % clipping boundary AND it is
    % not completely inside the
    % boundary - then its on the
    % boundary
    
    check_result = 1;
    
else
    check_result = 0;
    
end

end


function [x y z] = get_node_coords( nodes, node_num )

% Description
% ===========
%
% GET_NODE_COORDS() returns the X Y Z coordinates for a node based
% NOT on its row index, but rather it's actual node number as listed
% in the NODES array, column 1.
%
% History
% =======
%
% Dec 01 2010 -- Created -- MFH

node_numbers = nodes( : , 1 ); % node nums in column 1

node_row = find( node_numbers == node_num );

x = nodes( node_row, 2 ); % col 2 is x coord

y = nodes( node_row, 3 ); % col 3 is y coord

z = nodes( node_row, 4 ); % col 4 is z coord

end


function [ nodes, fibers ] = read_network( filename )

% Description
% ===========
%
% read_network() reads in a network file and returns a matrix of fibers
% and nodal coordinates in the network. Note: we need to make sure the
% network file does not have empty lines at the bottom of the file or we
% might have some problems.
%
% Output
% ======
%
% nodes is a N x 4 array of N nodes
% column 1 is the node number ordered least to greatest
% column 2 ix the x coord
% column 3 is the y coord
% column 4 is the z coord
%
% fibers is a M x 3 array of M fibers
% column 1 is the fiber number ordered as given in the vect file
% column 2 is the fiber node a
% column 3 is the fiber node b
%
% History
% =======
%
% Dec 01 2010 -- Modified from my del_gen_mesh() code from Oct 2010 -- MFH

raw_data = dlmread( filename, ' ', 1, 0); % load vect file skip header

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

