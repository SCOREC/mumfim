%% Examples on how to use network_generator
% Can create periodic delaunay, non periodic delaunay, periodic voronoi,
%   and non periodic voronoi
% Different networks may not work with very few points
% 100 seed points voronoi -> ~ 1000 nodes , 1500 elements
% 100 seed points delaunay -> ~ 600 nodes , 1000 elements

%% Create default network

% Periodic voronoi network with 50 seed points
network_generator('temp.txt');

%% Create periodic delaunay with 40 seed points and rve boundary at 0.6

nPoints = 40;
clip_boundary = 0.6;
network_generator('temp.txt',nPoints,clip_boundary,'delaunay');

%% Create non periodic voronoi network

nPoints = 30;
clip_boundary = 0.5;
networkType = 'voronoi';
periodic = false;

network_generator('temp.txt',nPoints,clip_boundary,networkType,periodic);

%% Create multiple random networks

n_networks = 3;
for i=1:n_networks
    network_generator(['temp_' num2str(i) '.txt']);
end
    
%% File input
% File is list of seed points in 3D
% Currently this is the only way to have non random or anisotropic networks
% If periodic then only supply points inside rve cube
% If non periodic then also supply points outside of rve cube

nPoints = 30; % Overwritten by number of points in file
clip_boundary = 0.5;
networkType = 'voronoi';
periodic = true;
fileInput = 'testInput.txt';

network_generator('temp.txt',nPoints,clip_boundary,networkType,periodic,fileInput);



