
% Case - Pure bending of beam, rollup. 
clear; 
clc; 
%% CONSTANTS. 
n_elems = 1; 
L0 = 1; 
dL = L0/n_elems;

x1 = [0, 0, 0]'; 
x2 = [L0, 0, 0]';
d3 = [1, 0, 0]'; % axial director. 
d2 = [0, 1, 0]'; % shear directors 2, 1. 
d1 = [0, 0, 1]'; 
dof = 12;  
fixed_dof = 12; % at first node. 
n_constraints = 6;
lambda = zeros(n_constraints,1);

lambdas = repmat(lambda, 1, n_elems+1); 
beam_dofs = repmat(dof, 1, n_elems +1);
beam_constraints = repmat(n_constraints, 1, n_elems + 1); 
beam_fixed_dofs = zeros(1, n_elems + 1);
beam_fixed_dofs(1) = fixed_dof;

[beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
beam_directors = repmat(vertcat(d1, d2, d3), 1, n_elems + 1);

beam_elems = create_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
beam = Beam(n_elems, beam_elems, dof, n_constraints); 

n_gauss_points = 1; % Gauss integration points. 

EA = 2; 
GA_Z = 2; 
GA_Y = 2; 
EI_Z = 2; 
EI_Y = 2; 
GJ = 2;

% material C matrix. 
C = diag([GA_Z, GA_Y, EA, EI_Z, EI_Y, GJ]);

% fixed and free nodes. 
fixed_nodes = [1]; 
free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, n_elems + 1);
%% Simulation. 
f_ext = zeros((beam.n_elems + 1) * beam.dof_per_node, 1); 

Fx = 0; 
Fy = 0; 
Fz = 0;

mx = 0;
my = 0; 
mz = 1;

% map spatial moments to directors. 


node2_x = (beam.dof_per_node )*( beam.n_nodes - 1) + 1; 
node2_y = node2_x + 1; 
node2_z = node2_x + 2;
node2_d1_mz = node2_z + 1; 
node2_d1_my = node2_z + 2; 
node2_d1_mx = node2_z + 3; 
node2_d2_mz = node2_z + 4; 
node2_d2_my = node2_z + 5; 
node2_d2_mx = node2_z + 6; 
node2_d3_mx = node2_z + 7; 
node2_d3_my = node2_z + 8; 
node2_d3_mz = node2_z + 9; 

% force at node position. 
f_ext(node2_x) = Fx;  
f_ext(node2_y) = Fy; 
f_ext(node2_z) = Fz; 

% moment applied.
f_ext(node2_d3_mz) = mz; 
f_ext(node2_d3_my) = my; 
f_ext(node2_d3_mx) = mx; 

% load steps. 
load_steps = 10; 

% tolerance. 
Tol = 1e-8; 
% number of gauss integration points. 
n_gauss_points = 1;

% iteration count. 
iter = 0; 
max_iter = 50; 

for i=1:load_steps
    force = (i/load_steps) * f_ext;
    [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, force, 1);
end

% plot beam. 
beam.show_config(L0 + 5, 5, 5, 1, "");
beam.display_end_node_pos();
beam.display_end_node_directors();