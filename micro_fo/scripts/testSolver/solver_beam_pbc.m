function solver_beam_pbc(filename,displacement)

%***********************************************************************
% Implements a network of Euler-Bernoulli 1D beam elements in a 3D cube
% This linear beam law accounts for axial, bending, and torsion forces.
%
% Periodic BCs
%
% Dan Fovargue
% April 2015
%***********************************************************************

% Get mesh
[nodes, fibers, clip_boundary, lrNodes, tbNodes, fbNodes] = read_periodic_network( filename );

% Initialize nodes and elements
nodes = nodes(:,2:4);
nodes = [nodes zeros(size(nodes,1),3)];
init_nodes = nodes;
dofspernode = 6;
elem = fibers(:,2:3);

plot_mesh(nodes,elem,1);

[leftb rightb topb bottomb frontb backb allb] = determineBoundaryNodes(nodes,clip_boundary);

% Hardcoded displacement of right RVE boundary
rvedisp = zeros(8,3);
rvedisp([2 4 6 8],1) = displacement;
nodes(:,1:3) = initialDisplace(nodes(:,1:3),rvedisp,clip_boundary);
disp_vector = reshape((nodes-init_nodes)',numel(nodes),1);
plot_mesh(nodes,elem,2);

% Material parameters
% E = 43200000;
% A = 3.8465*10^-8;
% G = 2500;

J = 2.51327*10^-7;
I1 = 1.25664*10^-7;
I2 = 1.25664*10^-7;

E = 6500;
G = 2500;
A = 0.00125663706;

ii=0;
mynorm = norm(disp_vector);

while(mynorm > 10^(-12) && ii<1)
    
    % Calculate stiffness matrix and forces
    lengths = calc_lengths(nodes,elem);
    stiffness_matrix = calc_matrix(nodes,elem,lengths,E,A,G,J,I1,I2);
    %     stiffness_matrix = stiffness_matrix + calc_geometric_matrix(nodes,elem,lengths,E,A,G,J,I1,I2,init_nodes);
    force_vector = -stiffness_matrix * disp_vector;
    
    % Apply BCs
    [stiffness_matrix force_vector] = apply_bcs(stiffness_matrix,force_vector, .... 
        leftb,rightb,topb,bottomb,frontb,backb,allb,lrNodes,tbNodes,fbNodes);
    
    % Solve
    disp_vector = stiffness_matrix\force_vector;
    
    mynorm = norm(disp_vector);
    fprintf('Norm = %d\n',mynorm);

    % Update nodes and displacements
    nodes = nodes + reshape(disp_vector,dofspernode,size(nodes,1))';
    disp_vector = reshape((nodes-init_nodes)',size(disp_vector,1),1);
    
%     ii = ii+1;
    
end

plot_mesh(nodes,elem,3);

% Recalculate force vector to have forces at boundaries
lengths = calc_lengths(nodes,elem);
stiffness_matrix = calc_matrix(nodes,elem,lengths,E,A,G,J,I1,I2);
force_vector = stiffness_matrix * disp_vector;

% Calculate modulus (using left and right faces - should be the same)
val = calcMod(rightb,leftb,force_vector);
fprintf('Modulus = %d  ,  %d\n',val(1)/displacement,val(2)/displacement);

end


function traction = calcMod(rightb,leftb,force_vector)

dpn = 6;
traction = [0 0];

tt1 = force_vector(rightb*dpn-5);
tt2 = force_vector(leftb*dpn-5);

traction(1) = sum(tt1);
traction(2) = sum(tt2);

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


function [stiffness_matrix force_vector] = ...
    apply_bcs(stiffness_matrix,force_vector,leftb,rightb,topb,bottomb,frontb,backb,allb,lrNodes,tbNodes,fbNodes)

dpn = 6;

% Fix position dof that corresponds to RVE side
for i=1:size(leftb,1)
    j = leftb(i)*dpn-5;
    force_vector(j) = 0.0;
    stiffness_matrix(j,:) = 0.0;
    stiffness_matrix(j,j) = 1.0;
end
for i=1:size(rightb,1)
    j = rightb(i)*dpn-5;
    force_vector(j) = 0.0;
    stiffness_matrix(j,:) = 0.0;
    stiffness_matrix(j,j) = 1.0;
end
for i=1:size(topb,1)
    j = topb(i)*dpn-4;
    force_vector(j) = 0.0;
    stiffness_matrix(j,:) = 0.0;
    stiffness_matrix(j,j) = 1.0;
end
for i=1:size(bottomb,1)
    j = bottomb(i)*dpn-4;
    force_vector(j) = 0.0;
    stiffness_matrix(j,:) = 0.0;
    stiffness_matrix(j,j) = 1.0;
end
for i=1:size(frontb,1)
    j = frontb(i)*dpn-3;
    force_vector(j) = 0.0;
    stiffness_matrix(j,:) = 0.0;
    stiffness_matrix(j,j) = 1.0;
end
for i=1:size(backb,1)
    j = backb(i)*dpn-3;
    force_vector(j) = 0.0;
    stiffness_matrix(j,:) = 0.0;
    stiffness_matrix(j,j) = 1.0;
end


% Set half of all boundary nodes' dofs to zero
% for u_i - u_j = 0 (periodic bc)
% First add the two equations together
for i=1:size(lrNodes,1)
    
    ldofs = [4 3 2 1 0];
    for d = ldofs
        j = lrNodes(i,1)*dpn-d;
        k = lrNodes(i,2)*dpn-d;
        force_vector(k) = force_vector(k) + force_vector(j);
        force_vector(j) = 0.0;
        stiffness_matrix(k,:) = stiffness_matrix(k,:) + stiffness_matrix(j,:);
        stiffness_matrix(j,:) = 0.0;
        stiffness_matrix(j,j) = 1.0;
        stiffness_matrix(j,k) = -1.0;
    end

end
for i=1:size(tbNodes,1)
    
    ldofs = [5 3 2 1 0];
    for d = ldofs
        j = tbNodes(i,1)*dpn-d;
        k = tbNodes(i,2)*dpn-d;
        force_vector(k) = force_vector(k) + force_vector(j);
        force_vector(j) = 0.0;
        stiffness_matrix(k,:) = stiffness_matrix(k,:) + stiffness_matrix(j,:);
        stiffness_matrix(j,:) = 0.0;
        stiffness_matrix(j,j) = 1.0;
        stiffness_matrix(j,k) = -1.0;
    end
    
end
for i=1:size(fbNodes,1)
    
    ldofs = [5 4 2 1 0];
    for d = ldofs
        j = fbNodes(i,1)*dpn-d;
        k = fbNodes(i,2)*dpn-d;
        force_vector(k) = force_vector(k) + force_vector(j);
        force_vector(j) = 0.0;
        stiffness_matrix(k,:) = stiffness_matrix(k,:) + stiffness_matrix(j,:);
        stiffness_matrix(j,:) = 0.0;
        stiffness_matrix(j,j) = 1.0;
        stiffness_matrix(j,k) = -1.0;
    end
    
end

end


function [stiffness_matrix] = calc_matrix(nodes,elem,lengths,E,A,G,J,I1,I2)

nnodes = size(nodes,1);
dpn = 6;
ndofs = nnodes*dpn;
stiffness_matrix = zeros(ndofs,ndofs);

for i=1:size(elem,1)
    
    L = lengths(i);
    n1 = elem(i,1);
    n2 = elem(i,2);
    
    cxx = (nodes(n2,1) - nodes(n1,1)) / lengths(i);
    cyx = (nodes(n2,2) - nodes(n1,2)) / lengths(i);
    czx = (nodes(n2,3) - nodes(n1,3)) / lengths(i);
    
    if cxx == 0 && cyx == 0
        if czx>0
            transform_sub_matrix = [0 0 1; 0 1 0; -1 0 0];
        else
            transform_sub_matrix = [0 0 -1; 0 1 0; 1 0 0];
        end
    else
        D = sqrt(cxx^2 + cyx^2);
        cxy = -cyx/D;
        cyy = cxx/D;
        czy = 0;
        cxz = -cxx*czx/D;
        cyz = -cyx*czx/D;
        czz = D;
    
        transform_sub_matrix = [cxx cyx czx; ...
                                cxy cyy czy; ...
                                cxz cyz czz ] ;
    end
    
    z = zeros(3,3);
    transform_matrix = [transform_sub_matrix z z z; ...
                        z transform_sub_matrix z z; ...
                        z z transform_sub_matrix z; ...
                        z z z transform_sub_matrix];
                    
    axial_matrix = zeros(12,12);
    axial_matrix(1,1) = 1.0;
    axial_matrix(1,7) = -1.0;
    axial_matrix(7,1) = -1.0;
    axial_matrix(7,7) = 1.0;
    axial_matrix = (E*A/L) * axial_matrix;
    
    torsion_matrix = zeros(12,12);
    torsion_matrix(4,4) = 1.0;
    torsion_matrix(4,10) = -1.0;
    torsion_matrix(10,4) = -1.0;
    torsion_matrix(10,10) = 1.0;
    torsion_matrix = (G*J/L) * torsion_matrix;
    
    bending_matrix = zeros(12,12);
    bending_matrix(2,2) = 12*E*I2/L^3;
    bending_matrix(3,3) = 12*E*I1/L^3;
    bending_matrix(5,5) = 4*E*I1/L;
    bending_matrix(6,6) = 4*E*I2/L;
    bending_matrix(8,8) = 12*E*I2/L^3;
    bending_matrix(9,9) = 12*E*I1/L^3;
    bending_matrix(11,11) = 4*E*I1/L;
    bending_matrix(12,12) = 4*E*I2/L;
    bending_matrix(2,6) = 6*E*I2/L^2;
    bending_matrix(3,5) = -6*E*I1/L^2;
    bending_matrix(5,9) = 6*E*I1/L^2;
    bending_matrix(6,8) = -6*E*I2/L^2;
    bending_matrix(8,12) = -6*E*I2/L^2;
    bending_matrix(9,11) = 6*E*I1/L^2;
    bending_matrix(2,8) = -12*E*I2/L^3;
    bending_matrix(3,9) = -12*E*I1/L^3;
    bending_matrix(5,11) = 2*E*I1/L;
    bending_matrix(6,12) = 2*E*I2/L;
    bending_matrix(2,12) = 6*E*I2/L^2;
    bending_matrix(3,11) = -6*E*I1/L^2;
    bending_matrix = bending_matrix + bending_matrix' - diag(diag(bending_matrix));
    
    k = axial_matrix + torsion_matrix + bending_matrix;
    
    local_matrix = transform_matrix' * k * transform_matrix;
    
    dofs_id = [n1*dpn-5 n1*dpn-4 n1*dpn-3 n1*dpn-2 n1*dpn-1 n1*dpn ...
               n2*dpn-5 n2*dpn-4 n2*dpn-3 n2*dpn-2 n2*dpn-1 n2*dpn];
    
    stiffness_matrix(dofs_id,dofs_id) = stiffness_matrix(dofs_id,dofs_id) + local_matrix;
    
end

end

function [stiffness_matrix] = calc_geometric_matrix(nodes,elem,lengths,E,A,G,J,I1,I2,init_nodes)

nnodes = size(nodes,1);
dpn = 6;
ndofs = nnodes*dpn;
stiffness_matrix = zeros(ndofs,ndofs);

for i=1:size(elem,1)
    
    L = lengths(i);
    n1 = elem(i,1);
    n2 = elem(i,2);
    
    cxx = (nodes(n2,1) - nodes(n1,1)) / lengths(i);
    cyx = (nodes(n2,2) - nodes(n1,2)) / lengths(i);
    czx = (nodes(n2,3) - nodes(n1,3)) / lengths(i);
    
    if cxx == 0 && cyx == 0
        if czx>0
            transform_sub_matrix = [0 0 1; 0 1 0; -1 0 0];
        else
            transform_sub_matrix = [0 0 -1; 0 1 0; 1 0 0];
        end
    else
        D = sqrt(cxx^2 + cyx^2);
        cxy = -cyx/D;
        cyy = cxx/D;
        czy = 0;
        cxz = -cxx*czx/D;
        cyz = -cyx*czx/D;
        czz = D;
    
        transform_sub_matrix = [cxx cyx czx; ...
                                cxy cyy czy; ...
                                cxz cyz czz ] ;
    end
    
    z = zeros(3,3);
    transform_matrix = [transform_sub_matrix z z z; ...
                        z transform_sub_matrix z z; ...
                        z z transform_sub_matrix z; ...
                        z z z transform_sub_matrix];
                    

    
    k = zeros(12,12);
    k(2,2) = 6/5;
    k(2,6) = L/10;
    k(2,8) = -6/5;
    k(2,12) = L/10;
    k(3,3) = 6/5;
    k(3,5) = -L/10;
    k(3,9) = -6/5;
    k(3,11) = -L/10;
    k(4,4) = J/A;
    k(4,10) = -J/A;
    k(5,5) = 2*L*L/15;
    k(5,9) = L/10;
    k(5,11) = -L*L/30;
    k(6,6) = 2*L*L/15;
    k(6,8) = -L/10;
    k(6,12) = -L*L/30;
    
    k(8,8) = 6/5;
    k(8,12) = -L/10;
    k(9,9) = 6/5;
    k(9,11) = L/10;
    k(10,10) = J/A;
    k(11,11) = 2*L*L/15;
    k(12,12) = 2*L*L/15;
    
    k = k + k' - diag(diag(k));
    d1 = (nodes(n2,1) - init_nodes(n2,1)) - (nodes(n1,1) - init_nodes(n1,1));
    d2 = (nodes(n2,2) - init_nodes(n2,2)) - (nodes(n1,2) - init_nodes(n1,2));
    d3 = (nodes(n2,3) - init_nodes(n2,3)) - (nodes(n1,3) - init_nodes(n1,3));
    k = (E*A/L * ( cxx*d1 + cyx*d2 + czx*d3 ))*k;
    
    local_matrix = transform_matrix' * k * transform_matrix;
    
    dofs_id = [n1*dpn-5 n1*dpn-4 n1*dpn-3 n1*dpn-2 n1*dpn-1 n1*dpn ...
               n2*dpn-5 n2*dpn-4 n2*dpn-3 n2*dpn-2 n2*dpn-1 n2*dpn];
    
    stiffness_matrix(dofs_id,dofs_id) = stiffness_matrix(dofs_id,dofs_id) + local_matrix;
    
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
clf;
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

function [ nodes, fibers, clip_boundary, lrNodes, tbNodes, fbNodes ] = ...
    read_periodic_network( filename )

% Copied from network clip function
% Except checks for optional clip boundary value included at end of file

nums = dlmread(filename, ' ', [0 0 0 5]);

raw_data = dlmread( filename, ' ', [1 0 nums(3) 8]); % load vect file skip header

lrNodes = dlmread( filename, ' ', [nums(3)+1 0 nums(3)+nums(4) 1]);
tbNodes = dlmread( filename, ' ', [nums(3)+nums(4)+1 0 nums(3)+nums(4)+nums(5) 1]);
fbNodes = dlmread( filename, ' ', [nums(3)+nums(4)+nums(5)+1 0 nums(3)+nums(4)+nums(5)+nums(6) 1]);

clip_boundary = dlmread( filename, ' ', [nums(3)+nums(4)+nums(5)+nums(6)+1 0 nums(3)+nums(4)+nums(5)+nums(6)+1 0]);

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

