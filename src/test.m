
clear;
clc; 
%% CONSTANTS. 
L = 1;
x1 = [0, 0, 0]'; 
x2 = [L, 0, 0]';
d3 = [1, 0, 0]'; % axial director. 
d2 = [0, 1, 0]'; % shear directors 2, 1. 
d1 = [0, 0, 1]'; 
dof = 12;  
n_constraints = 6;
lambda = zeros(n_constraints,1);
s = 0; 

fixed_dof = 1:12; % set first node fixed, for all kinematic variables.
fixed_const = 2*dof+1:2*dof+n_constraints; % fixed constraints at first node. 
all_dof = 1:2*dof + 2*n_constraints; % all dofs + constraints. 

% free degrees of freedom. 
free_dof = setdiff(all_dof, [fixed_dof, fixed_const]);
% material property matrix. 
C = eye(6);

% reference curvature and strains. 
gamma_ref = [0;0; 1];
omega_ref = [0;0;0];

% create beam class object. 
beam = Beam_elem(x1, x2, L, d3, d2, d1, d3, d2, d1, lambda, lambda, dof, n_constraints, fixed_dof, gamma_ref, omega_ref); 
%% Axial force. 
dof_per_node = 12; 
f_ext = zeros(2*dof_per_node, 1); 
Fx = 0; 
Fy = 2; 
Fz = 0;
Mx = 0;

% x force index. 
node2_x = dof + 1;

% tolerance. 
Tol = 1e-8; % 1e-10 
% number of gauss integration points. 
n_gauss_points = 1;

% iteration count. 
iter = 0; 
max_iter = 50; 

% delta_u_lambda = 1;
% main loop. 
verbose = 1; 

load_steps = 2; 
for i=1:load_steps
    F_x = (Fx/load_steps) * i;
    F_y = (Fy/load_steps) * i; 
    F_z = (Fz/load_steps) * i; 
    M_x = (Mx/load_steps) * i; 
    f_ext = zeros(2*dof_per_node, 1); 
    f_ext(node2_x) = F_x; 
    f_ext(node2_x + 1) = F_y;
    f_ext(node2_x + 2) = F_z; 
    f_ext(node2_x + 5) = 2* M_x;
    [iter] = Newtons_method(beam, n_gauss_points, C, gamma_ref, omega_ref, max_iter, Tol, free_dof, f_ext, n_constraints, verbose);
end

% beam.compute_f_int_u_displacement(u_tot, 1,C)

%% plot beam. 
x_lim = L + Fx + 1; 
y_lim = 4 + Fy + 1; 
z_lim = 4 + Fz + 1;
scaling = 0;
title = "Force at node 2: [Fx, Fy, Fz] = [" + num2str(F_x) + "N, " + num2str(F_y) + "N, " + num2str(F_z) + "N ]"; 
output1 = "Initial Node 2 position: [" + num2str(beam.x2(1)) + "m, " +  num2str(beam.x2(2)) + "m, " +  num2str(beam.x2(3)) + "m]";
output2 = "Node 2 position: [" + num2str(beam.x2_t(1)) + "m, " +  num2str(beam.x2_t(2)) + "m, " +  num2str(beam.x2_t(3)) + "m]";
beam.show_config(-x_lim, x_lim, -y_lim, y_lim, -z_lim, z_lim, scaling, title, 1);
disp(output1);
disp(output2);
beam.display_node_directors();

% check orhonormality condition. 
disp(transpose(beam.d1_B) * beam.d2_B); 
disp(transpose(beam.d1_B) * beam.d3_B);  
disp(transpose(beam.d2_B) * beam.d3_B); 



