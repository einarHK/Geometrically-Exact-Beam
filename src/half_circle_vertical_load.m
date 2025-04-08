
clear; 
clc; 
%% CONSTANTS. 
% initial beam length. 
L = 1; 
% number of elements. 
n_elems = 21; 
% kinematic degrees of freedom per node. 
dof = 12; 
% constraints per node. 
n_constraints = 6; 
% Initial lagrange multiplers for node. 
lambda = zeros(n_constraints, 1); 

% lambda multipliers for all beam elements. 
lambdas = repmat(lambda, 1, n_elems+1); 
% kinematic dof for each beam element. 
beam_dofs = repmat(dof, 1, n_elems +1);
% constraints for each beam element. 
beam_constraints = repmat(n_constraints, 1, n_elems + 1);
% the fixed dof for each beam element. 
beam_fixed_dofs = zeros(1, n_elems + 1);
% fix first and last node of the beam. 
fixed_dof = 12;
beam_fixed_dofs(1) = fixed_dof;
beam_fixed_dofs(n_elems + 1) = fixed_dof;  

% material property matrix. 
C = eye(6); 
C(1,1) = 75; 
C(2,2) = 75; 
C(3,3) = 100; 
C(4,4) = 100; 
C(5,5) = 100; 
C(6,6) = 200; 

% discretize the domain into beam coordinates, beam lengths and directors. 
% [beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
Radius = 2/pi;
bend_angle = pi/2; % 90 degree bend. 
[beam_coordinates, beam_directors, beam_lengths] = discretize_circle(Radius, bend_angle, n_elems);

% beam elements and beam object. 
beam_elems = create_curved_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
beam = Beam(n_elems, beam_elems, dof, n_constraints); 

% Gauss quadrature points. 
n_gauss_points = 1; 

% max iter count. 
max_iter = 50; 

% Tolerance. 
Tol = 1e-12; 

% dimension for e. 
e_dim = 6; 

% the fixed and free nodes. 
fixed_nodes = [1, n_elems + 1]; 
free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, n_elems + 1); 
% fixed_hq = beam.n_nodes * beam.dof_per_node + beam.constraint_per_node + 1: beam.n_nodes * beam.dof_per_node + 2 * beam.constraint_per_node;
% free_dof = setdiff(free_dof, fixed_hq); 
%% SIMULATION. 

% load steps. 
load_steps = 50; 

% load per node. 
fx = 0;
fy = 0; 
fz = 15; 

nodal_force_pos = [2, 3, 4, 9, 10, 11, 12, 13, 14, 15, 19, 20, 21];
f_ext = zeros(beam.n_nodes * beam.dof_per_node, 1); 

for i=2:beam.n_nodes - 1
    x_pos = (i-1) * beam.dof_per_node + 1; 
    y_pos = (i-1) * beam.dof_per_node + 2; 
    z_pos = (i-1) * beam.dof_per_node + 3; 
    f_ext(x_pos) = fx;
    f_ext(y_pos) = fy; 
    f_ext(z_pos) = fz;
end

f_ext = zeros(beam.n_nodes * beam.dof_per_node, 1);
for i=1:length(nodal_force_pos)
    node_pos = nodal_force_pos(i); 
    fx_pos = (node_pos - 1) * beam.dof_per_node + 1; 
    fy_pos = (node_pos - 1) * beam.dof_per_node + 2; 
    fz_pos = (node_pos - 1) * beam.dof_per_node + 3; 
    
    if node_pos < 8        
        f_ext(fx_pos) = -10; 
        f_ext(fz_pos) = -20; 
    elseif node_pos < 18
        f_ext(fx_pos) = 7.5; 
        f_ext(fy_pos) = -7.5; 
        f_ext(fz_pos) = 15; 
    else
        f_ext(fy_pos) = 10;
        f_ext(fz_pos) = -20;
    end

end

for i=1:load_steps
    force = f_ext * (i/load_steps); 
    % solve using Newton-Rhapson method. 
    [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, force, 1, free_dof); 
    
end


%% PLOTS. 

scale = 0.01; 
x_lim1 = -1; 
x_lim2 = 1; 
y_lim1 = -1; 
y_lim2 = 1; 
z_lim1 = -1; 
z_lim2 = 1;

beam.plot_undeformed_state(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, "");
hold on; 
beam.plot_deformed_beam(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, "Deformed (\color{red}red\color{black}) & Undeformed (\color{blue}blue\color{black})");
hold on;
coords = beam.get_deformed_beam_coords(); 
for i=1:length(nodal_force_pos)
    node_pos = nodal_force_pos(i); 
    coord = coords(:,node_pos); 
    scatter3(coord(1), coord(2), coord(3), color="yellow");
    hold on;
end

% plot node positions. 

for i=1:length(beam.beam_elements)
    beam_elem = beam.beam_elements(i);
    x_pos = beam_elem.x1_t(1); 
    y_pos = beam_elem.x1_t(2); 
    z_pos = beam_elem.x1_t(3); 
    dx = x_pos - beam_elem.x1(1); 
    dy = y_pos - beam_elem.x1(2); 
    dz = z_pos - beam_elem.x1(3); 
    d11 = beam_elem.d1_A(1); 
    d12 = beam_elem.d1_A(2); 
    d13 = beam_elem.d1_A(3); 
    d21 = beam_elem.d2_A(1); 
    d22 = beam_elem.d2_A(2); 
    d23 = beam_elem.d2_A(3);
    d31 = beam_elem.d3_A(1); 
    d32 = beam_elem.d3_A(2); 
    d33 = beam_elem.d3_A(3);
    output = "Node " + num2str(i) + " (x, y, z) = ( " + num2str(x_pos) +  ", " + num2str(y_pos) + ", " + num2str(z_pos) + " )"; 
    output2 = "Node " + num2str(i) + " (dx, dy, dz) = ( " + num2str(dx) +  ", " + num2str(dy) + ", " + num2str(dz) + " )"; 
    output3 = "Node Original Position: ( " + num2str(beam_elem.x1(1)) + ", " + num2str(beam_elem.x1(2)) + ", " + num2str(beam_elem.x1(3)) + " )";
    output4 = "Node Original director d1: ( " + num2str(d11) + ", " + num2str(d12) + ", " + num2str(d13) + " )"; 
    output5 = "Node Original director d1: ( " + num2str(d21) + ", " + num2str(d22) + ", " + num2str(d23) + " )"; 
    output6 = "Node Original director d1: ( " + num2str(d31) + ", " + num2str(d32) + ", " + num2str(d33) + " )"; 

    disp(output);
    disp(output2); 
    disp(output3); 
    disp(output4); 
    disp(output5);
    disp(output6);

    fprintf("\n"); 
end

x_pos = beam_elem.x2_t(1); 
y_pos = beam_elem.x2_t(2); 
z_pos = beam_elem.x2_t(3); 
dx = x_pos - beam_elem.x2(1); 
dy = y_pos - beam_elem.x2(2); 
dz = z_pos - beam_elem.x2(3); 
d11 = beam_elem.d1_B(1); 
d12 = beam_elem.d1_B(2); 
d13 = beam_elem.d1_B(3); 
d21 = beam_elem.d2_B(1); 
d22 = beam_elem.d2_B(2); 
d23 = beam_elem.d2_B(3);
d31 = beam_elem.d3_B(1); 
d32 = beam_elem.d3_B(2); 
d33 = beam_elem.d3_B(3);

output = "Node " + num2str(i+1) + " (x, y, z) = ( " + num2str(x_pos) +  ", " + num2str(y_pos) + ", " + num2str(z_pos) + " )"; 
output2 = "Node " + num2str(i+1) + " (dx, dy, dz) = ( " + num2str(dx) +  ", " + num2str(dy) + ", " + num2str(dz) + " )"; 
output3 = "Node Original Position: ( " + num2str(beam_elem.x2(1)) + ", " + num2str(beam_elem.x2(2)) + ", " + num2str(beam_elem.x2(3)) + " )";
output4 = "Node Original director d1: ( " + num2str(d11) + ", " + num2str(d12) + ", " + num2str(d13) + " )"; 
output5 = "Node Original director d2: ( " + num2str(d21) + ", " + num2str(d22) + ", " + num2str(d23) + " )"; 
output6 = "Node Original director d3: ( " + num2str(d31) + ", " + num2str(d32) + ", " + num2str(d33) + " )"; 
disp(output);
disp(output2); 
disp(output3); 
disp(output4); 
disp(output5); 
disp(output6); 

% calculate amount of asymmetry in z between nodes.
z_error = 0;
for i=1:length(beam.beam_elements)
    beam1 = beam.beam_elements(i); 
    i2 = beam.n_elements + 1 - i; 
    if i2 <= i
        break;
    end
    beam2 = beam.beam_elements(i2); 
    z1 = beam1.x1_t(3); 
    z2 = beam2.x2_t(3); 
    dz = z2 - z1;
    disp("Node " + num2str(i) + " and " + num2str(i2) + ": " + num2str(dz));
    z_error = z_error + abs(dz); 
    
end

z_error

