

% Quarter of a circular arc with total arc length of 1 m. 
clear; 
clc;
%% CONSTANTS. 
% initial beam length. 
L = 1; 
% number of elements. 
n_elems = 5; 
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

% directors. 
d1 = transpose([0, 0, 1]); 
d2 = transpose([0, 1, 0]); 
d3 = transpose([1, 0, 0]); 
beam_directors = repmat(vertcat(d1, d2, d3), 1, n_elems + 1);

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
bend_angle = (1/2) * pi; % 90 degree bend. 
x_start = [0, 0, 0]'; 
x_end = [L, 0, 0]';
% [beam_coordinates, beam_directors, beam_lengths] = discretize_circle(Radius, bend_angle, n_elems);
[beam_coordinates, beam_lengths] = discretize_domain(x_start, x_end, n_elems);

% beam elements and beam object. 
beam_elems = create_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
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
%% SIMULATION. 
load_steps = 5; 

% force vectors. 
fx = 0; 
fy = 0; 
fz = 50; 

% external force vector. 
f_ext = zeros((beam.n_nodes) * beam.dof_per_node, 1); 

for i=2:beam.n_nodes-1
    x_pos = (i-1) * beam.dof_per_node + 1;
    y_pos = (i-1) * beam.dof_per_node + 2; 
    z_pos = (i-1) * beam.dof_per_node + 3; 
    f_ext(x_pos) = fx; 
    f_ext(y_pos) = fy; 
    f_ext(z_pos) = fz;
end

beam.display_end_node_pos();
% iterate over each load step. 
for i=1:load_steps
    force = f_ext * (i/load_steps); 
    % solve using Newton-Rhapson method. 
    [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, force, 1, free_dof); 
end
%% PLOT RESULTS. 
scale = 0.01; 
x_lim1 = -0.1; 
x_lim2 = 1; 
y_lim1 = -1; 
y_lim2 = 1; 
z_lim1 = -5; 
z_lim2 = 5;

beam.plot_undeformed_state(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, "");
hold on; 
beam.plot_deformed_beam(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, "Deformed (\color{red}red\color{black}) & Undeformed (\color{blue}blue\color{black})");

% display nodal positions. 
for i=2:length(beam.beam_elements)
    beam_elem = beam.beam_elements(i); 
    x = beam_elem.x1_t(1); 
    y = beam_elem.x1_t(2); 
    z = beam_elem.x1_t(3); 
    output = "Node " + num2str(i) + ": (x, y, z) = " + "( " + num2str(x) + ", " + num2str(y) + ", " + num2str(z) + " )";
    disp(output);
end


