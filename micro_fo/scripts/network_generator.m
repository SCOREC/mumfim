function network_generator(outputFile,varargin)
%********************************************************
% INPUT
% =====
%  Mandatory inputs
%  ----------------
%    outputFile - Name of filename to output network into.
%
%  Optional inputs
%  ---------------
%    numPoints - Number of points to seed the network.
%                If periodic this is exact, if not then a multiple of this
%                is used to seed a larger cube surrounding the desired
%                region.
%    clip_boundary - Cube boundary value, will be from -clip_boundary to
%                    +clip_boundary in each direction.
%    networkType - 'delaunay' or 'voronoi'
%    periodic - true or false for periodic network
%    inputFile - Name of file to get seed points from. No checking is done
%                on these points so its up to you to make sure they make
%                sense. Should be ascii file with N rows and 3 columns.
%                Periodicity should still work because it just makes copies
%                of the data in the original cube.
%
% DESCRIPTION
% ===========
% Generates a delaunay or voronoi based fiber network.
% Can be periodic or not.
% Uses randomly distributed seed points unless input file is specified.
% For uniform or anisotropic networks use input file with such points.
%
% This does contain code to plot the network but its skipped.
%
% HISTORY
% =======
%
% 4/2015 - Big overhaul to clean up and add support for periodic networks
%          Changed input, output, structure, comments, etc
%          - Dan Fovargue
%
% 9/2014 - Rewritten to support delaunay and voronoi network - Dan Fovargue
%
%
%
% 2/21 - Began re-writing find_edges--> temporary code starting line 92.
%        once completed, this will be saved in a separate (likely multiple
%        separate) function(s). 
%       
% 3/4  - Problem from lines 215-249. Not saving the correct segments on each
%        face. Added to the SVN group repository.
%
% 3/7  - Fixed problem from lines 215-249. 
%
% 3/8  - find_edges now a function called in the main calling section
%
% 3/17 - Begun breaking donw find_edges into smaller functions:
%        find_facets, find_planar_facets, connected_tris .
%
% 3/21 - Completed find_facets which outputs a structure in which each
%        element contains the nodes of the triangle facets on a given 
%        voronoi cell. 
%
% 3/23 - Began constructing find_planar_facets which should create a
%        structure that contains each of the nodes of the triangles that 
%        make up a given facet in one element.
%
% 4/1  - Spawned a new branch to ONLY extract the fiber network from
%        a Voronoi tesselation
%********************************************************

% Set default values
numPoints = 50;
clip_boundary = 0.5;
networkType = 'voronoi';
periodic = true;
fromFile = false;
inputFile = '';

plot_flag = false;

% Check and assign inputs
if nargin<1
    error('Input error: Must specify output file');
elseif nargin>7
    error('Input error: Too many inputs');
end
if nargin>=2
    numPoints = varargin{1};
end
if nargin>=3
    clip_boundary = varargin{2};
end
if nargin>=4
    networkType = varargin{3};
    if ~(strcmp(networkType,'delaunay') || strcmp(networkType,'voronoi'))
        error('Input error: Network type must be string, either ''delaunay'' or ''voronoi''');
    end
end
if nargin>=5
    periodic = varargin{4};
end
if nargin>=6
    inputFile = varargin{5};
    fromFile = true;
end

tic

fprintf('Generating network with the following parameters:\n');
fprintf('outputFile = %s\n',outputFile);
fprintf('numPoints = %i\n',numPoints);
fprintf('clip_boundary = %f\n',clip_boundary);
fprintf('networkType = %s\n',networkType);
fprintf('periodic = %i\n',periodic);
fprintf('fromFile = %i\n',fromFile);
fprintf('inputFile = %s\n',inputFile);

% Make cube_length depending on periodic-ness
if periodic
    % For periodic we only need random points inside the final cube
    cube_length = 2*clip_boundary;
else
    % If not periodic make cube larger than final cube so we get fibers
    % attached to the boundary
    nonPeriodicLengthMultiplier = 2;
    cube_length = nonPeriodicLengthMultiplier*2*clip_boundary;
    numPoints = numPoints * nonPeriodicLengthMultiplier^3;
end

% Points to seed the network
if fromFile
    
    pts = load( inputFile );
    numPoints = size(pts,1);
    
    if size(pts,2) ~= 3
        error('Data from input file is the wrong size');
    end
    
    x = pts(:,1);
    y = pts(:,2);
    z = pts(:,3);
    
else
    
    volume_v = cube_length*rand(numPoints,3) - cube_length/2.0;
    
    x = volume_v(:,1)';
    y = volume_v(:,2)';
    z = volume_v(:,3)';
    
end

% If periodic then create copies of points around the originals
if periodic
    periodicOffsets = [ 1  0  0;  0  1  0;  0  0  1;
                       -1  0  0;  0 -1  0;  0  0 -1;
                        1  1  0;  1 -1  0; -1  1  0;
                       -1 -1  0;  1  0  1;  1  0 -1;
                       -1  0  1; -1  0 -1;  0  1  1;
                        0  1 -1;  0 -1  1;  0 -1 -1; 
                        1  1  1;  1  1 -1;  1 -1  1; 
                        1 -1 -1; -1  1  1; -1  1 -1; 
                       -1 -1  1; -1 -1 -1];
                   
    periodicOffsets = periodicOffsets*2*clip_boundary;
    
    x0 = x;
    y0 = y;
    z0 = z;
    
    for i=1:size(periodicOffsets,1)
        x = [x; x0+periodicOffsets(i,1)];
        y = [y; y0+periodicOffsets(i,2)];
        z = [z; z0+periodicOffsets(i,3)];
    end
end

% Create initial network
if strcmp(networkType,'delaunay')
    nodes = [x(:),y(:),z(:)];
    [cells] = delaunayn( nodes ); % Generates delaunay cells
    cells = mat2cell(cells, ones(1,size(cells,1)) );
else
    [nodes,cells] = voronoin( [x(:),y(:),z(:)] ); % Generates voronoi points
end

% By deafult this is skipped
if plot_flag
    plotInitNetwork(x,y,z); 
end

if strcmp(networkType,'voronoi')
    % Remove cells that are connected to node 1 the infinite node
    % (Not applicable for delaunay network)
    cells = filter_cells( cells );
end

% Generate point cloud for each voronoi cell
sharp_edges = [];

for n = 1 : length( cells )
    
    unique_vertices_in_bubble = cells{n};
    node_xyz = nodes( unique_vertices_in_bubble, : );
    
    % CONVHULLN returns the triangular facets of the cube faces
    node_vertices = convhulln( node_xyz );
    
    % TRIREP is used to find the FEATUREDGES of the cube faces
    node_tri = TriRep( node_vertices, node_xyz );
    local_edges = featureEdges( node_tri, pi/1000 ); % pi/1000 is our bend threshold
    
    % CONVERT LOCAL TO GLOBAL INDEXING
    for i = 1 : size( local_edges, 1 )
        global_edges(i,1) =  unique_vertices_in_bubble( local_edges(i,1) );
        global_edges(i,2) =  unique_vertices_in_bubble( local_edges(i,2) );
    end
    
    sharp_edges = [ sharp_edges;
                    global_edges ]; 
                
end

% Gather unique fibers
fibers = unique( sharp_edges, 'rows' );

% Some fibers get doubled up, at least with delaunay networks
% Get rid of these
fibers = removeDoubledFibers(fibers);

% Renumber nodes
[nodes, fibers, cells] = renumber_nodes( nodes, fibers, cells );

% Put IDs in front of nodes and fibers
nodes = [(1:size(nodes,1))' nodes];
fibers = [(1:size(fibers,1))' fibers];

% Clip network to boundary
[nodes, fibers] = network_n_clip_translate(nodes,fibers,clip_boundary);

% Finalize the network
if periodic
    
    % Check and label periodicity
    [leftb rightb topb bottomb frontb backb] = determineBoundaryNodes(nodes,clip_boundary);
    l = clip_boundary*2;
    eps = 10^-4;
    lrNodes = checkPeriodicBoundaryNodes(rightb,leftb,[l 0 0],nodes,eps);
    tbNodes = checkPeriodicBoundaryNodes(topb,bottomb,[0 l 0],nodes,eps);
    fbNodes = checkPeriodicBoundaryNodes(frontb,backb,[0 0 l],nodes,eps);
    
    % Output along with periodic connection information
    write_mesh_periodic(outputFile,nodes,fibers,clip_boundary,lrNodes,tbNodes,fbNodes);
    
else

    % Remove any solo fibers
    [nodes, fibers] = removeIndependentFibers(nodes,fibers);
    
    % Output network
    write_mesh(outputFile,nodes,fibers,clip_boundary);
    
end

fprintf('Finished generating network. Total time: %f\n',toc);

% Plot final network - skipped by default
if plot_flag
    plot_mesh(nodes(:,2:4),fibers(:,2:3),2);
end

end
% End of main network_generator function


% Begin support functions...


function fibers = removeDoubledFibers(fibers)

if size(fibers,2)~=2
    error('removeDoubledFibers expects 2 columns in fibers')
end

numFibers = size(fibers,1);

% Make sure node a < node b
for i=1:numFibers
    if fibers(i,1) > fibers(i,2)
        temp = fibers(i,1);
        fibers(i,1) = fibers(i,2);
        fibers(i,2) = temp;
    elseif fibers(i,1) == fibers(i,2)
        error('Huge error - fiber from node a to node a')
    end
end

fibers = unique(fibers,'rows');

end


function [nnodes] = checkPeriodicBoundaryNodes(b1,b2,dx,nodes,eps)

% find i,j where:  b1(i) - dx = b2(j)

if length(b1) ~= length(b2)
    error('Number of nodes on opposite boundaries aren''t equal')
end

nnodes = zeros(length(b1),2);
nn = 1;

for i=1:length(b1)
    
    ii = b1(i);
    p1 = nodes(ii,2:4);
    
    for j=1:length(b2)
        
        jj = b2(j);
        p2 = nodes(jj,2:4);
        
        % Check if two points are "equal" i.e. separated by dx
        if norm((p1-dx)-p2,1)<eps
            nnodes(nn,:) = [ii jj];
            nn = nn + 1;
            break;
        end
    end
end

% Make sure nothing weird happened
if ~( isequal(b1,nnodes(:,1)) && isequal(b2,sort(nnodes(:,2))) )
    error('Boundaries don''t match up')
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


function [leftb rightb topb bottomb frontb backb] ...
    = determineBoundaryNodes(nodes,clip_boundary)

leftb = [];
rightb = [];
topb = [];
bottomb = [];
frontb = [];
backb = [];

nodes = nodes(:,2:4);

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

end


function write_mesh_periodic( filename, nodes, fibers, clip_boundary, bn_lr, bn_tb, bn_fb)

total_nodes         = size( nodes, 1 );     % count rows = num nodes
total_fibers        = size( fibers, 1 );    % count rows = num fibers

tot_deg_freedom     = 3 * total_nodes;      % x y z values for nodes

fileid = fopen(filename,'w'); % writes over any existing file of this name

fprintf(fileid,'%i %i %i %i %i %i\n', ...
    total_nodes, tot_deg_freedom, total_fibers, size(bn_lr,1), size(bn_tb,1), size(bn_fb,1));

% Output data for each fiber
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

% Output periodic boundary relations
bn = [bn_lr; bn_tb; bn_fb];
for i = 1:length(bn)
    fprintf(fileid,'%i %i\n',bn(i,1),bn(i,2));
end

% Output cube boundary
fprintf(fileid,'%f\n',clip_boundary);

fclose( fileid );

end


function write_mesh( filename, nodes, fibers, clip_boundary)

total_nodes         = size( nodes, 1 );     % count rows = num nodes
total_fibers        = size( fibers, 1 );    % count rows = num fibers

tot_deg_freedom     = 3 * total_nodes;      % x y z values for nodes

fileid = fopen(filename,'w'); % writes over any existing file of this name

fprintf(fileid,'%i %i %i\n', total_nodes, tot_deg_freedom, total_fibers);

% Output data for each fiber
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

% Output cube boundary
fprintf(fileid,'%f\n',clip_boundary);

fclose( fileid );

end


function [new_cells] = filter_cells( cells )

% INPUT
% =====
%
% OUTPUT
% ======
% NEW_CELLS - all cells except ones with node 1 numbered 1...N cells
%
% HISTORY
% =======
% Mar 2011 -- Last mod

new_cell_count = 1;

for n = 1 : length( cells )  
   
    if cells{n}(1) ~= 1; % Avoid the cells that have the infinite node number 1
        
        new_cells{ new_cell_count } = cells{ n };
        
        new_cell_count = new_cell_count + 1;
       
    end
    
end

end


function [nodes,fibers,cells] = renumber_nodes( nodes, fibers, cells )

% New version of this does whatever the old one did expect it doesn't
% actually output the network to file. 4/2015 Dan Fovargue

% DESCRIPTION
% ===========
% Function writes out a network where some instances in NODES do not
% appear in the FIBERS array -- such as in the case of a Voronoi network
% where we skip cells linked to an infinite node -- renumbering the nodes
% that actually do appear in the network.
%
% INPUT
% =====
% NODES - m-by-3 array of xyz coordinates for m nodes indexed 1-m
% FIBERS - m-by-2 array of start-end node indices for 1-m fibers
%
% OUTPUT
% ======
% NODES - resized to remove unused nodes
% FIBERS - nodes re-numbered sequentially
%
% HISTORY
% =======
% Apr 2011 -- Modified from write_network() -- MFH
% Jun 2011 -- Modified to update indexing in "cells" -- JDH

total_fibers = size( fibers, 1 ); % num rows = num fibers regardless of indexing

% Collect old node indices
old_node_indices = unique( fibers(1:end) ); % 1 x N array -- sorted least-greatest
total_nodes = length( old_node_indices ); % number of nodes in play

% Renumber FIBERS node indices and collect new NODES array
for n = 1 : total_fibers
    
    old_start_node       = fibers(n,1);
    
    old_end_node         = fibers(n,2);
    
    new_start_index      = find( old_node_indices == old_start_node );
    
    new_end_index        = find( old_node_indices == old_end_node );
    
    fibers(n,1)          = new_start_index;
    
    fibers(n,2)          = new_end_index;
    
end

nodes = nodes( old_node_indices, : );

total_cells = length(cells);

for n = 1 : total_cells
    
    total_nodes_in_cell = length(cells{n});
    
    for m = 1 : total_nodes_in_cell
        cells{n}(m) = find( old_node_indices == cells{n}(m) );
    end
    
end

end


% This function is outdated but its only display anyways
function [] = plotInitNetwork(x,y,z)

% Plot of all the seed points
figure('Color',[1 1 1]);

scatter3( x, y, z, '.', 'r' );
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])

% Plot of only the points in the -0.5 to 0.5 range
j = 1;
for i=1:length(x)
    if x(i) >= -0.5 && x(i) <= 0.5
        if y(i) >= -0.5 && y(i) <= 0.5
            if z(i) >= -0.5 && z(i) <= 0.5
                points(j,1:3) = [x(i); y(i); z(i)];
                j = j+1;
            end
        end
    end
end

figure('Color',[1 1 1]);
plot3(points(:,1),points(:,2),points(:,3),'MarkerSize',20,'Marker','.','LineStyle','none');
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
grid off
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'ZTickLabel',[])
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);

% Put a box around fibers
Xfmin = -0.5;
Xfmax = 0.5;

Yfmin = -0.5;
Yfmax = 0.5;

Zfmin = -0.5;
Zfmax = 0.5;

% yplusface
[X,Z]=meshgrid(linspace(Xfmin,Xfmax,2),linspace(Zfmax,Zfmin,2));
Y=Yfmax*ones(2);
hand_yplus=surface(X,Y,Z);
set(hand_yplus,'FaceAlpha',0.1);
set(hand_yplus,'FaceColor','none','LineWidth',2,'EdgeColor',[1 0 0]);

% yminusface
[X,Z]=meshgrid(linspace(Xfmin,Xfmax,2),linspace(Zfmax,Zfmin,2));
Y=Yfmin*ones(2);
hand_yminus=surface(X,Y,Z);
set(hand_yminus,'FaceAlpha',0.1);
set(hand_yminus,'FaceColor','none','LineWidth',2,'EdgeColor',[1 0 0]);

% xplusface
[Y,Z]=meshgrid(linspace(Yfmin,Yfmax,2),linspace(Zfmax,Zfmin,2));
X=Xfmax*ones(2);
hand_xplus=surface(X,Y,Z);
set(hand_xplus,'FaceAlpha',0.1);
set(hand_xplus,'FaceColor','none','LineWidth',2,'EdgeColor',[1 0 0]);

% xminusface
[Y,Z]=meshgrid(linspace(Yfmin,Yfmax,2),linspace(Zfmax,Zfmin,2));
X=Xfmin*ones(2);
hand_xminus=surface(X,Y,Z);
set(hand_xminus,'FaceAlpha',0.1);
set(hand_xminus,'FaceColor','none','LineWidth',2,'EdgeColor',[1 0 0]);

% Zplusface
[X,Y]=meshgrid(linspace(Xfmin,Xfmax,2),linspace(Yfmin,Yfmax,2));
Z=Zfmax*ones(2);
hand_zplus=surface(X,Y,Z);
set(hand_zplus,'FaceAlpha',0.1);
set(hand_zplus,'FaceColor','none','LineWidth',2,'EdgeColor',[1 0 0]);

% Zminusface
[X,Y]=meshgrid(linspace(Xfmin,Xfmax,2),linspace(Yfmin,Yfmax,2));
Z=Zfmin*ones(2);
hand_zminus=surface(X,Y,Z);
set(hand_zminus,'FaceAlpha',0.1);
set(hand_zminus,'FaceColor','none','LineWidth',2,'EdgeColor',[1 0 0]);

axis image
view(130,15);

csvwrite('points_out.csv',points);

end

