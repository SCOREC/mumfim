function [nodes, disp_vector] = solver_truss_phantom(filename,displacement,deform_type)

%***********************************************************************
% Implements a network of trusses in a 3D cube.
% Truss network is made up of a main network and a phantom network. The
% motivation is to remove the singularity in the stiffness matrix that
% arises for Voronoi networks by adding a "phantom" Delaunay network that
% is much more compliant.
%
% Victor Chan
% Aug 2015
%***********************************************************************

% Get mesh
[nodes, fibers, clip_boundary] = read_network( filename );

% Initialize nodes and elements
nodes = nodes(:,2:4);
nnodes = size(nodes,1);
elem = fibers(:,2:3);

plot_mesh(nodes,elem,1);

[leftb rightb topb bottomb frontb backb allb] = determineBoundaryNodes(nodes,clip_boundary);

% Material parameters
E = 43200000;
radius = 3.49911e-08;
A = pi*(radius)^2;
B = 1.2; %constant for nonlinear fiber force calculation.

% Vector of 1s and 0s - 1s correspond to fixed DOFs
fixedDOFs = getFixedDOFs(nodes,leftb, rightb, topb, bottomb, frontb, backb, allb);

% Hardcoded displacement of right RVE boundary
rvedisp = zeros(8,3);
init_lengths = calc_lengths(nodes,elem);

rvedisp([2 4 6 8],1) = displacement;
% rvedisp([1 3 5 7],1) = 0.0;
nodes = initialDisplace(nodes,rvedisp,clip_boundary);

lengths = calc_lengths(nodes,elem);
forces = calc_forces(elem,lengths,init_lengths,E,A,B,deform_type);
force_vector = calc_force_vector(nodes,elem,lengths,forces);
force_vector = force_bcs(force_vector,fixedDOFs);

mynorm = norm(force_vector);
fprintf('Norm = %d\n',mynorm);
    
if strcmp(deform_type,'linear') == 1
    stiffness_matrix = calc_matrix(nodes,elem,lengths,E,A); 
    stiffness_matrix = matrix_bcs(stiffness_matrix,fixedDOFs);
end

if strcmp(deform_type,'nonlinear') == 1
    stiffness_matrix = calc_nonlin_matrix(nodes,elem,lengths,init_lengths,forces,E,A,B);
    stiffness_matrix = matrix_bcs(stiffness_matrix,fixedDOFs);
end

plot_mesh(nodes,elem,2);
ii=0;

while(mynorm > 10^(-12) && ii<20)

    disp_vector = -stiffness_matrix\force_vector; 

    nodes = nodes + reshape(disp_vector,3,nnodes)'; %equivalent to update_nodes() in biotissue code.

    lengths = calc_lengths(nodes,elem);

    forces = calc_forces(elem,lengths,init_lengths,E,A,B,deform_type);

    force_vector = calc_force_vector(nodes,elem,lengths,forces);

    force_vector = force_bcs(force_vector,fixedDOFs);

    mynorm = norm(force_vector);

    fprintf('Norm = %d\n',mynorm);

    if strcmp(deform_type,'linear') == 1
        stiffness_matrix = calc_matrix(nodes,elem,lengths,E,A); 
        stiffness_matrix = matrix_bcs(stiffness_matrix,fixedDOFs);
    end

    if strcmp(deform_type,'nonlinear') == 1
        stiffness_matrix = calc_nonlin_matrix(nodes,elem,lengths,init_lengths,forces,E,A,B);
        stiffness_matrix = matrix_bcs(stiffness_matrix,fixedDOFs);
    end

    ii=ii+1;
end

plot_mesh(nodes,elem,3);

% Recalculate force vector to have forces at boundaries
lengths = calc_lengths(nodes,elem);
forces = calc_forces(elem,lengths,init_lengths,E,A,B,deform_type);
force_vector = calc_force_vector(nodes,elem,lengths,forces);

% sum of forces on face (using left and right faces - should be the same)
val = calcMod(rightb,leftb,force_vector);
fprintf('xfxp = %e \n',val(1));
fprintf('xfxn = %e \n',val(2));

% Calculate microscale stresses
stress = calc_stress(nodes,force_vector,allb);
fprintf('microscale stresses = %e, %e, %e, %e, %e, %e\n', ...
        stress(1), stress(2), stress(3), stress(4), stress(5), stress(6));
end

%% Functions used in main program. 
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

function fixedDOFs = getFixedDOFs(nodes,leftb, rightb, topb, bottomb, frontb, backb, allb)

dpn = 3;
fixedDOFs = zeros(size(nodes,1)*dpn,1);

% All boundary nodes fixed in x, y, z
    for i=1:size(allb,1)

        fixedDOFs(allb(i)*dpn-2) = 1;
        fixedDOFs(allb(i)*dpn-1) = 1;
        fixedDOFs(allb(i)*dpn) = 1;

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

function [forces] = calc_forces(elem,lengths,init_lengths,E,A,B,deform_type)

    forces = zeros(size(elem,1),1);
    if strcmp(deform_type, 'linear') == 1
        for i=1:size(elem,1)
            forces(i) = E*A*(lengths(i)/init_lengths(i) - 1);
        end
    end

    if strcmp(deform_type, 'nonlinear') == 1
        for i=1:size(elem,1)
            epsf = 0.5*( (lengths(i)/init_lengths(i))^2-1.0 ); % Equivalent to green_strain in biotissue code.
            forces(i) = (E*A/B)*(exp(B*epsf)-1.0);
        end
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

        force_vector(n1*3 - 2) = force_vector(n1*3 - 2) - cx*f;
        force_vector(n1*3 - 1) = force_vector(n1*3 - 1) - cy*f;
        force_vector(n1*3 - 0) = force_vector(n1*3 - 0) - cz*f;
        
        force_vector(n2*3 - 2) = force_vector(n2*3 - 2) + cx*f;
        force_vector(n2*3 - 1) = force_vector(n2*3 - 1) + cy*f;
        force_vector(n2*3 - 0) = force_vector(n2*3 - 0) + cz*f;

    end

end

function [force_vector] = force_bcs(force_vector,fixedDOFs)

for i=1:length(force_vector)
    if(fixedDOFs(i))
        force_vector(i) = 0.0;
    end
end

end

function [stiffness_matrix] = calc_matrix(nodes,elem,lengths,E,A)

nnodes = size(nodes,1);
ndofs = nnodes*3;
stiffness_matrix = zeros(ndofs,ndofs);

    for i=1:size(elem,1)

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

end

function [stiffness_matrix] = calc_nonlin_matrix(nodes,elem,lengths,init_lengths,forces,E,A,B)

nnodes = size(nodes,1);
ndofs = nnodes*3;
stiffness_matrix = zeros(ndofs,ndofs);

    for i=1:size(elem,1)

        f = forces(i);
        l = lengths(i);
        l0 = init_lengths(i);
        
        n1 = elem(i,1);
        n2 = elem(i,2);

        cx = (nodes(n2,1) - nodes(n1,1)) / lengths(i);
        cy = (nodes(n2,2) - nodes(n1,2)) / lengths(i);
        cz = (nodes(n2,3) - nodes(n1,3)) / lengths(i);

        epsf = 0.5*( (l/l0)^2-1.0 ); % equivalent to green strain in biotissue code.
        k = (A*E*l/l0^2)*exp(B*epsf)-f/l;  % equivalent to dfdl-fl in biotissue code.
%         k = 0.5*epsf*l/l0^2*f-f/l;

        transform_sub_matrix = [k*cx^2+f/l  k*cx*cy k*cx*cz; ...
                                k*cx*cy k*cy^2+f/l  k*cy*cz; ...
                                k*cx*cz k*cy*cz k*cz^2+f/l ] ;
        local_matrix = [ transform_sub_matrix -transform_sub_matrix; ...
                        -transform_sub_matrix  transform_sub_matrix ];

        dofs_id = [n1*3-2 n1*3-1 n1*3 n2*3-2 n2*3-1 n2*3];

        for r=1:6
            for c=1:6
                stiffness_matrix(dofs_id(r),dofs_id(c)) = ...
                    stiffness_matrix(dofs_id(r),dofs_id(c)) + local_matrix(r,c);
            end
        end


    end

end

function [stiffness_matrix] = matrix_bcs(stiffness_matrix,fixedDOFs)

for i=1:length(fixedDOFs)
    if(fixedDOFs(i))
        
        stiffness_matrix(i,:) = 0.0;
        stiffness_matrix(i,i) = 1.0;
        
    end
end

end

function plot_mesh( nodes, fibers, fign)

% INPUT
% =====
% NODES     N-by-3 array of xyz coordinates for N nodes
% FIBERS    N-by-2 array of start-end nodes for N fibers

figure(fign);
clf;
plot3( nodes(:,1) , nodes(:,2), nodes(:,3),'o', ...  
    'LineWidth',1, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor', [0.5 0.5 0.8], 'MarkerSize',10 );

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
xlabel('X'); ylabel('Y'); zlabel('Z');
% view([0,0])
axis equal;

end

%% Post solver functions.
function [stress] = calc_stress(nodes,force_vector,allb)
% [allb*6-5, allb*6-4, allb*6-3]
% [force_vector(allb*6-5), force_vector(allb*6-4), force_vector(allb*6-3)]
% [nodes(allb,1), nodes(allb,2), nodes(allb,3)]
    dpn = 3;
    stress = zeros(1,6);
    xfx = sum(nodes(allb,1).*force_vector(allb*dpn - 2));
    yfx = sum(nodes(allb,2).*force_vector(allb*dpn - 2));
    zfx = sum(nodes(allb,3).*force_vector(allb*dpn - 2));
    
    xfy = sum(nodes(allb,1).*force_vector(allb*dpn - 1));
    yfy = sum(nodes(allb,2).*force_vector(allb*dpn - 1));
    zfy = sum(nodes(allb,3).*force_vector(allb*dpn - 1));
    
    xfz = sum(nodes(allb,1).*force_vector(allb*dpn));
    yfz = sum(nodes(allb,2).*force_vector(allb*dpn));
    zfz = sum(nodes(allb,3).*force_vector(allb*dpn));
    
    stress(1) = xfx;
    stress(2) = 0.5*(xfy + yfx);
    stress(3) = 0.5*(zfx + xfz);
    stress(4) = yfy;
    stress(5) = 0.5*(yfz + zfy);
    stress(6) = zfz;
   
end

function [stress] = calc_stress_sum_face(nodes,force_vector,allb)
% [allb*6-5, allb*6-4, allb*6-3]
% [force_vector(allb*6-5), force_vector(allb*6-4), force_vector(allb*6-3)]
% [nodes(allb,1), nodes(allb,2), nodes(allb,3)]
    stress = zeros(1,6);
    xfx = sum(nodes(allb,1).*force_vector(allb*6 - 5));
    yfx = sum(nodes(allb,2).*force_vector(allb*6 - 5));
    zfx = sum(nodes(allb,3).*force_vector(allb*6 - 5));
    
    xfy = sum(nodes(allb,1).*force_vector(allb*6 - 4));
    yfy = sum(nodes(allb,2).*force_vector(allb*6 - 4));
    zfy = sum(nodes(allb,3).*force_vector(allb*6 - 4));
    
    xfz = sum(nodes(allb,1).*force_vector(allb*6 - 3));
    yfz = sum(nodes(allb,2).*force_vector(allb*6 - 3));
    zfz = sum(nodes(allb,3).*force_vector(allb*6 - 3));
    
    stress(1) = xfx;
    stress(2) = 0.5*(xfy + yfx);
    stress(3) = 0.5*(zfx + xfz);
    stress(4) = yfy;
    stress(5) = 0.5*(yfz + zfy);
    stress(6) = zfz;
   
end

function traction = calcMod(rightb,leftb,force_vector)

dpn = 3;
traction = [0 0];

tt1 = force_vector(rightb*dpn-2);
tt2 = force_vector(leftb*dpn-2);

traction(1) = sum(tt1);
traction(2) = sum(tt2);

end