% 45 degree bend - cantilever beam case. 
clear;
clc;
%% CONSTANTS. 
n_elems = 20; % number of beam elements. 

dof = 12; % kinematic degrees of freedom per node. 
fixed_dof = 12; % fix the kinematic dof at the first node. 
n_constraints = 6; % number of constraints per node. 
lambda = zeros(n_constraints,1); % lambda multipliers per node - initialized as zeros. 

% lambda multipliers for all beam elements. 
lambdas = repmat(lambda, 1, n_elems+1); 
% kinematic dof for each beam element. 
beam_dofs = repmat(dof, 1, n_elems +1);
% constraints for each beam element. 
beam_constraints = repmat(n_constraints, 1, n_elems + 1);
% the fixed dof for each beam element. 
beam_fixed_dofs = zeros(1, n_elems + 1);
beam_fixed_dofs(1) = fixed_dof;

% discretize the domain into beam coordinates, beam lengths and directors. 
% [beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
Radius = 100;
bend_angle = (45/180) * pi; % 45 degree bend. 
[beam_coordinates, beam_directors, beam_lengths] = discretize_circle(Radius, bend_angle, n_elems);

% beam elements and beam object. 
beam_elems = create_curved_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
beam = Beam(n_elems, beam_elems, dof, n_constraints); 

% Material properties. 
A = 1; 
I = (1/12); 
J = (1/6);
E = 1e7; 

% Material property matrix. 
GJ = E*I;
EI2 = E*I; 
EI3 = E*I;
EA = E*A; 
GA2 = 5e6; 
GA3 = 5e6; 
C = diag([GA2, GA3, EA, EI2, EI3, GJ]);  

% Gauss quadrature points. 
n_gauss_points = 1; 

% max iter count. 
max_iter = 50; 

% fixed and free nodes. 
fixed_nodes = [1]; 
free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, n_elems + 1); 

% Tolerance. 
Tol = 1e-6; 
%% SIMULATION. 
load_steps = 10; 
% force components - end node. 
fx = 0; 
fy = 0; 
fz = 600; 
% force node positions. 
node_fx = beam.n_elements * beam.dof_per_node + 1; 
node_fy = beam.n_elements * beam.dof_per_node + 2;
node_fz = beam.n_elements * beam.dof_per_node + 3; 

% iterate over each load step. 
for i=1:load_steps
    f_ext = zeros((beam.n_nodes) * beam.dof_per_node, 1); 
    f_ext(node_fx) = fx * (i/load_steps); 
    f_ext(node_fy) = fy * (i/load_steps); 
    f_ext(node_fz) = fz * (i/load_steps);
    % solve using Newton-Rhapson method. 
    [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, f_ext, 1, free_dof); 
end
%% PLOT RESULTS. 
x_lim1 = -50; 
x_lim2 = 50; 
y_lim1 = 0; 
y_lim2 = 50; 
z_lim1 = -20; 
z_lim2 = 60; 

scaling = 0; 
scale_factor = 2;
title = ""; 
beam.plot(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scaling, title, scale_factor);

fprintf(1, '\n');

disp("Numerical solution: ")
node_indices = [n_elems + 1]; 
beam_coords = beam.get_deformed_beam_coords(); 
end_node_coord = beam_coords(:, n_elems + 1); 
x1 = end_node_coord(1); 
y1 = end_node_coord(2); 
z1 = end_node_coord(3);
disp("Load level: " + num2str(fz) + " - End node coords: (" + num2str(x1) + ", " + num2str(y1) + ", " + num2str(z1) + ")")

fprintf(1, '\n');
% analytical solution - Cardona. 
disp("Solution from Cardona: ")
display_45deg_analytical_solution();
