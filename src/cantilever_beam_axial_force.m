
% Cantilever beam - Axial force case. 
clear;
clc; 
%% CONSTANTS. 
n_elems = 10; % number of beam elements. 
L = 1; % initial length. 
x1 = transpose([0, 0, 0]); % start coordinate. 
x2 = transpose([L, 0, 0]); % end coordinate. 

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
d1 = transpose([0, 0, 1]); 
d2 = transpose([0, 1, 0]); 
d3 = transpose([1, 0, 0]); 
d_vector = [d1;d2;d3]; 
beam_directors = repmat(d_vector, n_elems + 1); 

% discretize the domain into beam coordinates, beam lengths and directors. 
[beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
beam_elems = create_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);

% beam elements and beam object. 
beam = Beam(n_elems, beam_elems, dof, n_constraints); 
%% SIMULATION. 
% max tolerance. 
TOL = 1e-12; 

% number of load steps. 
load_steps = 10; 

% force components. 
fx = 1; 
fy = 0; 
fz = 0; 

% material property matrix. 
C = diag([1, 1, 1, 1, 1, 1]); 

% max number of iterations. 
max_iter = 50; 

% gauss integration points. 
n_gauss_points = 1; 

% fixed and free nodes. 
fixed_nodes = [1]; 
free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, n_elems + 1); 

for i=1:load_steps
    f_ext = zeros((beam.n_elements + 1) * beam.dof_per_node, 1); 
    fx_node = beam.n_elements * beam.dof_per_node + 1; 
    fy_node = beam.n_elements * beam.dof_per_node + 2; 
    fz_node = beam.n_elements * beam.dof_per_node + 3; 
    f_ext(fx_node) = (i/load_steps) * fx; 
    f_ext(fy_node) = (i/load_steps) * fy; 
    f_ext(fz_node) = (i/load_steps) * fz; 
    [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, TOL, f_ext, 1, free_dof); 

end

%% PLOTS. 
x_lim1 = 0; 
x_lim2 = fx + L; 
y_lim1 = -1; 
y_lim2 = 1 + fy; 
z_lim1 = 0; 
z_lim2 = 2 + fz;
scale = 0.1; 

title_str =  "Axial Elongation: Fx = " + num2str(fx) +  "N " + ", Fy = " + num2str(fy) + "N " + ", Fz = " + num2str(fz) + "N"; 
beam.plot_deformed_beam(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, title_str);


