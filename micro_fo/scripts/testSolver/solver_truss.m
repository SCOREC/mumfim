function [stiffness] = solver_truss(filename,displacement)

%***********************************************************************
% Implements a network of 1D truss elements in a 3D cube
%
% (As expected, will work for delaunay network but not for voronoi)
%
% (Was set up to mimic the nonlinear solve embedded in the mutliscale model
% so its more complicated than a linear truss solve needs to be)
%
% Dan Fovargue
% Oct 2014
%***********************************************************************

[nodes, fibers, clip_boundary] = read_network( filename );

% Initialize nodes and elements
nodes = nodes(:,2:4);
nnodes = size(nodes,1);
ndofs = nnodes*3;
elem = fibers(:,2:3);
nelem = size(elem,1);

rvedisp = zeros(8,3);

plot_mesh(nodes,elem,1);

[leftb rightb topb bottomb frontb backb allb] = determineBoundaryNodes(nodes,clip_boundary);

fixedDOFs = getFixedDOFs(nodes,leftb, rightb, topb, bottomb, frontb, backb, allb);

init_lengths = calc_lengths(nodes,elem);
init_nodes = nodes;

% displacement = clip_boundary*2*(0.1);
rvedisp([2 4 6 8],1) = displacement;
% rvedisp(2,1) = displacement*10;
nodes = initialDisplace(nodes,rvedisp,clip_boundary);

plot_mesh(nodes,elem,2);

E = 6500;
A = 0.001256637;

lengths = calc_lengths(nodes,elem);

forces = calc_forces(elem,lengths,init_lengths,E,A);

force_vector = calc_force_vector(nodes,elem,lengths,forces);

force_vector = force_bcs(force_vector,fixedDOFs);

mynorm = norm(force_vector);

fprintf('Norm = %d\n',mynorm);

stiffness_matrix = calc_matrix(nodes,elem,lengths,forces,E,A,fixedDOFs);

ii=0;
while(mynorm > 10^(-8) && ii<1)
    
    disp_vector = stiffness_matrix\force_vector;
    
    nodes = nodes + reshape(disp_vector,3,nnodes)';
    
    lengths = calc_lengths(nodes,elem);
    
    forces = calc_forces(elem,lengths,init_lengths,E,A);
    
    force_vector = calc_force_vector(nodes,elem,lengths,forces);
    
    force_vector = force_bcs(force_vector,fixedDOFs);
    
    mynorm = norm(force_vector);
    
    fprintf('Norm = %d\n',mynorm);
    
    stiffness_matrix = calc_matrix(nodes,elem,lengths,forces,E,A,fixedDOFs);
    
    ii=ii+1;
end

% Find stiffness
lengths = calc_lengths(nodes,elem);
forces = calc_forces(elem,lengths,init_lengths,E,A);
force_vector = calc_force_vector(nodes,elem,lengths,forces);
xp_traction = traction(force_vector,rightb,1);
stiffness = -xp_traction/displacement;

plot_mesh(nodes,elem,3);

end


function total = traction(force_vector,bound,dir)

total = 0.0;
for i=1:size(bound,1)
    total = total + force_vector(bound(i)*3 - (3-dir));
end

end


function nodes = initialDisplace(nodes,rvedisp,clip_boundary)

for i=1:size(nodes,1)
    xx = (nodes(i,1)+clip_boundary)/(2*clip_boundary);
    yy = (nodes(i,2)+clip_boundary)/(2*clip_boundary);
    zz = (nodes(i,3)+clip_boundary)/(2*clip_boundary);
    
    xd = ( ((rvedisp(1,1))*(1-xx)+(rvedisp(2,1))*(xx))*(1-yy) + ...
           ((rvedisp(5,1))*(1-xx)+(rvedisp(6,1))*(xx))*(yy) )*(1-zz) + ...
         ( ((rvedisp(3,1))*(1-xx)+(rvedisp(4,1))*(xx))*(1-yy) + ...
           ((rvedisp(7,1))*(1-xx)+(rvedisp(8,1))*(xx))*(yy) )*(zz);
       
    yd = ( ((rvedisp(1,2))*(1-xx)+(rvedisp(2,2))*(xx))*(1-yy) + ...
           ((rvedisp(5,2))*(1-xx)+(rvedisp(6,2))*(xx))*(yy) )*(1-zz) + ...
         ( ((rvedisp(3,2))*(1-xx)+(rvedisp(4,2))*(xx))*(1-yy) + ...
           ((rvedisp(7,2))*(1-xx)+(rvedisp(8,2))*(xx))*(yy) )*(zz);
       
    zd = ( ((rvedisp(1,3))*(1-xx)+(rvedisp(2,3))*(xx))*(1-yy) + ...
           ((rvedisp(5,3))*(1-xx)+(rvedisp(6,3))*(xx))*(yy) )*(1-zz) + ...
         ( ((rvedisp(3,3))*(1-xx)+(rvedisp(4,3))*(xx))*(1-yy) + ...
           ((rvedisp(7,3))*(1-xx)+(rvedisp(8,3))*(xx))*(yy) )*(zz);
       
    nodes(i,1) = nodes(i,1) + xd;
    nodes(i,2) = nodes(i,2) + yd;
    nodes(i,3) = nodes(i,3) + zd;
end

end

function [lengths] = calc_lengths(nodes,elem)

lengths = zeros(size(elem,1),1);
for i=1:size(elem,1)
    lengths(i) = sqrt( (nodes(elem(i,2),1)-nodes(elem(i,1),1))^2 + ...
                       (nodes(elem(i,2),2)-nodes(elem(i,1),2))^2 + ... 
                       (nodes(elem(i,2),3)-nodes(elem(i,1),3))^2 );

end

end

function fixedDOFs = getFixedDOFs(nodes,leftb, rightb, topb, bottomb, frontb, backb, allb)

dpn = 3;
fixedDOFs = zeros(size(nodes,1)*dpn,1);

% All boundary nodes fixed in x, y, z
for i=1:size(allb,1)
    
    fixedDOFs(allb(i)*dpn-2) = 1;
    fixedDOFs(allb(i)*dpn-1) = 1;
    fixedDOFs(allb(i)*dpn  ) = 1;
    
end

end

function [forces] = calc_forces(elem,lengths,init_lengths,E,A)

forces = zeros(size(elem,1),1);
for i=1:size(elem,1)
    forces(i) = E*A*(lengths(i)/init_lengths(i) - 1);
end

end


function [force_vector] = calc_force_vector(nodes,elem,lengths,forces)

nnodes = size(nodes,1);
ndofs = nnodes*3;
force_vector = zeros(ndofs,1);

for i=1:size(elem,1)
    
    f = forces(i);
    n1 = elem(i,1);
    n2 = elem(i,2);
    
    cx = (nodes(n2,1) - nodes(n1,1)) / lengths(i);
    cy = (nodes(n2,2) - nodes(n1,2)) / lengths(i);
    cz = (nodes(n2,3) - nodes(n1,3)) / lengths(i);
    
    force_vector(n2*3 - 2) = force_vector(n2*3 - 2) - cx*f;
    force_vector(n2*3 - 1) = force_vector(n2*3 - 1) - cy*f;
    force_vector(n2*3 - 0) = force_vector(n2*3 - 0) - cz*f;
    
    force_vector(n1*3 - 2) = force_vector(n1*3 - 2) + cx*f;
    force_vector(n1*3 - 1) = force_vector(n1*3 - 1) + cy*f;
    force_vector(n1*3 - 0) = force_vector(n1*3 - 0) + cz*f;
    
end

end

function [force_vector] = force_bcs(force_vector,fixedDOFs)

for i=1:length(force_vector)
    if(fixedDOFs(i))
        force_vector(i) = 0.0;
    end
end

end

function [stiffness_matrix] = calc_matrix(nodes,elem,lengths,forces,E,A,fixedDOFs)

nnodes = size(nodes,1);
ndofs = nnodes*3;
stiffness_matrix = zeros(ndofs,ndofs);

for i=1:size(elem,1)
    
    f = forces(i);
    n1 = elem(i,1);
    n2 = elem(i,2);
    
    cx = (nodes(n2,1) - nodes(n1,1)) / lengths(i);
    cy = (nodes(n2,2) - nodes(n1,2)) / lengths(i);
    cz = (nodes(n2,3) - nodes(n1,3)) / lengths(i);
    
    transform_sub_matrix = [cx^2  cx*cy cx*cz; ...
                            cx*cy cy^2  cy*cz; ...
                            cx*cz cy*cz cz^2 ] ;
    local_matrix = E*A/lengths(i) * ...
                   [ transform_sub_matrix -transform_sub_matrix; ...
                    -transform_sub_matrix  transform_sub_matrix ];
    
    dofs_id = [n1*3-2 n1*3-1 n1*3 n2*3-2 n2*3-1 n2*3];
                
    for r=1:6
        for c=1:6
            stiffness_matrix(dofs_id(r),dofs_id(c)) = ...
                stiffness_matrix(dofs_id(r),dofs_id(c)) + local_matrix(r,c);
        end
    end

    
end

for i=1:length(fixedDOFs)
    if(fixedDOFs(i))
        
        stiffness_matrix(i,:) = 0.0;
        stiffness_matrix(i,i) = 1.0;
        
    end
end

end



function [leftb rightb topb bottomb frontb backb allb] ...
    = determineBoundaryNodes(nodes,clip_boundary)

leftb = [];
rightb = [];
topb = [];
bottomb = [];
frontb = [];
backb = [];

for i=1:size(nodes,1)
    
    if(nodes(i,1)==clip_boundary)
        rightb = [rightb; i];
    elseif(nodes(i,1)==-clip_boundary)
        leftb = [leftb; i];
    end
    
    if(nodes(i,2)==clip_boundary)
        topb = [topb; i];
    elseif(nodes(i,2)==-clip_boundary)
        bottomb = [bottomb; i];
    end
    
    if(nodes(i,3)==clip_boundary)
        frontb = [frontb; i];
    elseif(nodes(i,3)==-clip_boundary)
        backb = [backb; i];
    end
    
end

allb = unique(sort([rightb; leftb; topb; bottomb; frontb; backb]));

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

