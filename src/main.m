%% INIT. 
clear; 
clc; 
%% CONSTANTS
L = 1;
x1 = [0, 0, 0]'; 
x2 = [L, 0, 0]';
d1 = [1, 0, 0]'; 
d2 = [0, 1, 0]';
d3 = [0, 0, 1]';
dof = 12;  

C = eye(6); % C material properties matrix,
n_constraints = 6; % number of constraints. 
lambda = zeros(n_constraints,1); % lagrange coefficients. 
fixed_dof = 1:12; % fixed degrees of freedom. 
beam = Beam(x1, x2, L, d1, d2, d3, d1, d2, d3, dof, lambda, lambda, fixed_dof, n_constraints); % beam element - 2 nodes. 

% reference strains. 
gamma_ref = [L;0;0]; 
omega_ref = [0;0;0]; 

s_vector = [0, 0, 0, 0, 0, 0]';
max_iter = 100; 
TOL_q = 1e-10;

q = beam.compute_q(0);
disp(q);

%% Test. 
s = 0;
H = beam.compute_H(s); 

disp(size(H));

%% Example - Beam Elongation. 
Fx = 1; 
f_ext = zeros(2 * beam.dof, 1); 
node2_x = beam.dof + 1; 
f_ext(node2_x) = Fx; 

% interpolation value s. 
s = 0;

% gauss quadrature integration points. 
n_gauss_points = 1; 

% plot_beam(beam);

% disp(size(Kt));
[iter] = Newtons_method(beam, 0, C, gamma_ref, omega_ref, max_iter, TOL_q, TOL_q, f_ext, n_gauss_points); 

% e = beam.compute_e(q, gamma_ref, omega_ref);
% B = beam.compute_B(s); 
% f = B' * C * e; 
% disp(f);

disp(["Total iterations: ", num2str(iter)]);
disp(["Force Fx = ", num2str(Fx), "N"]);
beam.display_node_coords();
plot_beam(beam);
