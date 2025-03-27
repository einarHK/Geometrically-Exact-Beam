
clear;
clc;

%% CONSTANTS. 
n_elems = 2; 
L0 = 1; 

x1 = [0, 0, 0]'; 
x2 = [L0, 0, 0]';
d3 = [1, 0, 0]'; % axial director. 
d2 = [0, 1, 0]'; % shear directors 2, 1. 
d1 = [0, 0, 1]'; 
dof = 12;  
fixed_dof = 12; % at first node. 
n_constraints = 6;
lambda = zeros(n_constraints,1);

gamma_ref = [0; 0; 1]; 
omega_ref = [0; 0; 0];

% lambda vector for each node. 
lambdas = repmat(lambda, 1, n_elems+1); 
beam_dofs = repmat(dof, 1, n_elems +1);
beam_constraints = repmat(n_constraints, 1, n_elems + 1); 
beam_fixed_dofs = zeros(1, n_elems + 1);
beam_fixed_dofs(1) = fixed_dof;

% calculate beam coordinates and beam lengths. 
[beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
beam_directors = repmat(vertcat(d1, d2, d3), 1, n_elems + 1);

% create beam elements and beam class object. 
beam_elems = create_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
beam = Beam(n_elems, beam_elems, dof, n_constraints);

beam1 = Beam_elem(x1, x2, L0, d3, d2, d1, d3, d2, d1, lambda, lambda, dof, n_constraints, fixed_dof, gamma_ref, omega_ref); 

C = eye(6);
% S1 = beam1.compute_S_mat(1, C, gamma_ref, omega_ref); 
% S2 = beam.compute_S_mat(1, C); 


