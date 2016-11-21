% DESCRIPTION
% ===========
% Script to animate network step iterations in the microscale problem.
%
% HISTORY
% =======
% Apr 04 2011 -- Created -- MFH
% Jul 14 2011 -- Updated to plot from raw network file -- MFH




function single_network_plotter_printer(fname,NUM_NETWORKS,CLIP_BOUNDARY)
       
    
for n = 1:NUM_NETWORKS


    % LOAD NETWORK
    % ============
    
    %fname = 'FCL_211rcg';
    input_file      = sprintf('%s_%i.txt',fname,n);
    output_file     = sprintf('%s_%i.jpg',fname,n);
         
    [nodes,fibers]      = read_network( input_file );
    
    nodes_array_1       = nodes;
    
    fibers_array_1      = fibers;
    
    num_nodes_array_1   = size( nodes, 1 ); % count rows to get num nodes

    
    
    
    % FIBER DATA (COLOR) ARRAY
    % ========================
    
    color_val = 1.5;
    
    data_array_1 = color_val * ones( 1, size(fibers_array_1,1) ) ; % rand or ones ok here
    
    simple_network_3D_plotter( nodes_array_1, fibers_array_1, data_array_1 , CLIP_BOUNDARY)
    
   
        
    
    % IF YOU WOULD LIKE TO PRINT THE PLOT TO A JPG
    % ============================================
    
    if 1 == 0
                
        
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 3]); 
        print('-djpeg100',output_file,'-r200');
        
        
        
    end
    
end




end








function [node_xyz, fiber_nodes] = read_network( filename )


% OUTPUT
% ======
% NODE_XYZ      N-by-3 array of xyz coordinates for N nodes
% FIBER_NODES   N-by-2 array of start-end nodes for N fibers


raw_data = dlmread( filename, ' ', 1, 0); % load vect file skip header

fiber_number                = raw_data( : , 1 );
fiber_node_a                = raw_data( : , 2 );
fiber_node_b                = raw_data( : , 3 );

fiber_nodes                 = [fiber_node_a fiber_node_b]; % <--- nodal connectivity

fiber_total                 = length( fiber_number );

fiber_node_a_x              = raw_data( : , 4 );
fiber_node_a_y              = raw_data( : , 5 );
fiber_node_a_z              = raw_data( : , 6 );

fiber_node_b_x              = raw_data( : , 7 );
fiber_node_b_y              = raw_data( : , 8 );
fiber_node_b_z              = raw_data( : , 9 );                      

node_stack                  = [ fiber_node_a ...      % node number
                                fiber_node_a_x ...    % node x coord                 
                                fiber_node_a_y ...    % node y coord                     
                                fiber_node_a_z;       % node z coord
                                
                                fiber_node_b ...
                                fiber_node_b_x ...
                                fiber_node_b_y ...
                                fiber_node_b_z ];
                            
node_stack                  = unique( node_stack , 'rows' );

node_xyz                    = node_stack(: , 2:4); % <--- nodal coordinates

end








function plot_mesh( nodes, fibers, step, iter )


% INPUT
% =====
% NODES     N-by-3 array of xyz coordinates for N nodes
% FIBERS    N-by-2 array of start-end nodes for N fibers


figure;



plot3( nodes(:,1) , nodes(:,2), nodes(:,3),'o', ...  
    'LineWidth',1, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor', [0.5 0.5 0.8], 'MarkerSize',3 );

%axis_val = max( nodes(1:end) ) + 0.1;     % find a max/min for plot

%axis( [-axis_val axis_val -axis_val axis_val -axis_val axis_val] );

%axis( [-1 +1 -1 +1 -1 +1] );

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

%axis( [-3 +3 -3 +3 -3 +3] );

title_val = sprintf('Step %i Iter %i',step, iter);

title(title_val,'FontSize',8);

end





function [] = simple_network_3D_plotter( nodes, fibers, stretch, CLIP_BOUNDARY )


% Description
% ===========
%
% Plots out a fiber network from a set of rendered 3D cylinders.
%
% Some of this code was taken form a post by Michael Garrity at Mathworks.
%
% INPUT
% =====
% NODES - N-by-3 array of xyz coordinates for N nodes
% FIBERS - N-by-2 array of start-end nodes for N fibers
% STRETCH - N-by-1 array of N fiber stretchs
%
% History
% =======
%
% Jan 2008 -- Created -- EAS 
% Dec 2010 -- Modified -- MFH


plot_fibers_switch = 1;
plot_nodes_switch = 0;
plot_seed_points_switch = 0;
fiber_box_switch = 1;


figure;

NumFib          = size( fibers,1 ); % rows = number of fibers

%Ax = Ax';
%Ay = Ay';
%Az = Az';

%pt1_vec = [Ax(:,1) Ay(:,1) Az(:,1)];
%pt2_vec = [Ax(:,2) Ay(:,2) Az(:,2)];

%fib_vec = pt2_vec-pt1_vec;

[x,y,z]         = cylinder(0.01,80); % ORIGINAL = (0.01, 80)


if plot_fibers_switch
    
for i = 1:NumFib
    
    
        node_A_index = fibers( i, 1 ); % col 1 = node A
        node_B_index = fibers( i, 2 ); % col 2 = node B
        
        node_A_xyz = nodes( node_A_index, : ); % row vector of x y z for node A
        node_B_xyz = nodes( node_B_index, : ); % row vector of x y z for node B
        
        fib_vector = node_B_xyz - node_A_xyz; % row vector for the fiber
        
        fib_len     = norm( fib_vector );
        
        fib_vector = fib_vector / fib_len; % normalize
        
        tf = hgtransform;
        
        hold on;

        
        
        %fib_vector  = fib_vec(i,:)/fib_len;

        % Cylinder is oriented along the z-axis,i.e. [0 0 1].
        % Use cross product to find rotation angle. Vectors as columns
        
        rot_axis=cross([0 0 1],fib_vector);
        
        % Check that the fiber is not already aligned along z
        if norm(rot_axis) > eps
            
            rot_angle = asin(norm(rot_axis));
            
            if (dot([0 0 1],fib_vector) < 0)
                rot_angle = pi-rot_angle;
            end
            
        else
            
            rot_axis = [0 0 1];
            rot_angle = 0;
        
        end
        
        %fc = 100 * rand(1,1) * ones(size(x)); % use only 1 random color for each fiber based on a random value 0-64
        
        
        %fc = 100 * stretch(i) * ones(size(x)); % Color fiber via fiber stretch
        
        fc = stretch(i) * ones(size(x));
        
        
        %fc = fcolor(i)*ones(size(x));

        surface(x,y,z,fc,'parent',tf);

        set(tf,'matrix',makehgtform('translate',node_A_xyz,...
            'axisrotate',rot_axis,rot_angle,...
            'scale',[1 1 fib_len]));
    
    
 
end
end



if plot_nodes_switch
    
    plot3(  nodes(:,1), ...                   % plot nodes
            nodes(:,2), ...
            nodes(:,3), ...
            '.', ...
            'LineWidth', 6, ...
            'MarkerEdgeColor','k', ...
            'MarkerFaceColor', [0 0 0], ...
            'MarkerSize', 20 );
    
end


if plot_seed_points_switch
    
    points = csvread('points_out.csv');
    plot3(points(:,1),points(:,2),points(:,3),'MarkerSize',20,'Marker','.','LineStyle','none');
end


hold off

%view(135,20);

view(130,15);

shading interp;
camlight;
lighting gouraud;

colormap( 'gray' ); % jet winter bone



%colorbar; % <--------------------- COLORBAR



set(gcf, 'color', 'white');
axis equal; 




%axis( [-3 +3 -3 +3 -3 +3] ); % 3D axis settings for a stretch in X
%axis( [-1 +1 -1 +1 -1 +1] ); 
axis( [-CLIP_BOUNDARY CLIP_BOUNDARY -CLIP_BOUNDARY CLIP_BOUNDARY -CLIP_BOUNDARY CLIP_BOUNDARY] ); 

label_axis_switch = 1;

if label_axis_switch

    xlabel('X','FontSize',24);
    ylabel('Y','FontSize',24);
    zlabel('Z','FontSize',24);

end


caxis( [0.0 2.0] );


set(gca,'FontSize',10,'LineWidth',1);

set(gca,'Visible','on')



if fiber_box_switch

    % Put a box around fibers
    
    Xfmin = min( nodes( : , 1) );
    Xfmax = max( nodes( : , 1) );
    
    Yfmin = min( nodes( : , 2) );
    Yfmax = max( nodes( : , 2) );
    
    Zfmin = min( nodes( : , 3) );
    Zfmax = max( nodes( : , 3) );
    
    
    
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
    

end



end