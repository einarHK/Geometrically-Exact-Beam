
clear;
clc; 

%% testing
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

free_dof = setdiff(all_dof, [fixed_dof, fixed_const]);
C = eye(6);

gamma_ref = [0;0; L];
omega_ref = [0;0;0];

beam = Beam_elem(x1, x2, L, d3, d2, d1, d3, d2, d1, lambda, lambda, dof, n_constraints, fixed_dof, gamma_ref, omega_ref); 
% f_int = beam.compute_f_int(1, C, gamma_ref, omega_ref); 
%% Axial force. 
dof_per_node = 12; 
f_ext = zeros(2*dof_per_node, 1); 
Fx = 2; 
Fy = 0; 
Fz = 0;

node2_x = dof + 1; 
f_ext(node2_x) = Fx;  

% tolerance. 
Tol = 1e-6; 
% number of gauss integration points. 
n_gauss_points = 1;

% residual = beam.compute_residual(n_gauss_points, C, gamma_ref, omega_ref, f_ext); 
% rhs = beam.compute_rhs(n_gauss_points, C, gamma_ref, omega_ref, f_ext);
s_vector = beam.compute_s(s, C, gamma_ref, omega_ref); 
S = beam.compute_S_mat(n_gauss_points, C, gamma_ref, omega_ref);

% iteration count. 
iter = 0; 
max_iter = 50; 

% free degrees of freedom. 
free_dof_res = beam.compute_free_dof();
free_dof_h = n_constraints+1:2*n_constraints; 
dof_h_q = 2 * beam.n_constraints;

fy = Expand_vector(free_dof, dof_h_q); 
fx = Expand_vector(free_dof, dof_h_q); 

% compute initial residual. 
residual = beam.compute_residual(n_gauss_points, C, gamma_ref, omega_ref, f_ext); 
residual = residual(free_dof_res); 

% compute constraint vector. 
h_q = beam.compute_h_q(); 
h_q = h_q(free_dof_h); 

% compute errors. 
error_residual = MSE_error(residual); 
error_h_q = MSE_error(h_q); 

% main loop. 
while (((error_residual > Tol) && (error_h_q > Tol)) || (iter < max_iter))
    delta_u = zeros(2*beam.dof, 1); % solution vector - u values. 
    delta_lambda = zeros(2*beam.n_constraints, 1); % solution vector - lambda values. 
    delta_u_lambda = [delta_u; delta_lambda]; 

    % compute S matrix. 
    % s_vector = beam.compute_s(s, C, gamma_ref, omega_ref); 
    S = beam.compute_S_mat(n_gauss_points, C, gamma_ref, omega_ref); 
    S = S(free_dof, free_dof); 

    % compute force terms. 
    f_int = beam.compute_f_int(n_gauss_points, C, gamma_ref, omega_ref); 
    f_H1 = beam.compute_f_H1(); 
    f_H2 = beam.compute_f_H2();
    f_H = [f_H1; f_H2]; 
    residual = f_int + f_H - f_ext; 
    % residual = residual(free_dof); 
    h_q = beam.compute_h_q();
    % solve linearized system. 
    rhs = [-residual; -h_q];
    rhs = rhs(free_dof); 
    delta_u_lambda(free_dof) = S \ rhs; 
    
    % disp(S);

    % update values. 
    i1 = 1; 
    i2 = 2 * beam.dof; 
    i3 = length(delta_u_lambda); 
    dx1 = delta_u_lambda(1:3); 
    delta_d1_A = delta_u_lambda(4:6); 
    delta_d2_A = delta_u_lambda(7:9); 
    delta_d3_A = delta_u_lambda(10:12); 
    dx2 = delta_u_lambda(13:15); 
    delta_d1_B = delta_u_lambda(16:18); 
    delta_d2_B = delta_u_lambda(19:21); 
    delta_d3_B = delta_u_lambda(22:24); 
    delta_u = delta_u_lambda(i1:i2); 
    delta_lambda = delta_u_lambda(i2+1:i3);
    delta_lambda_A = delta_lambda(1:beam.n_constraints); 
    delta_lambda_B = delta_lambda(beam.n_constraints + 1: 2*beam.n_constraints);

    % update beam values. 
    beam.update_params(dx1, dx2, delta_d1_A, delta_d2_A, delta_d3_A, delta_d1_B, delta_d2_B, delta_d3_B); 
    beam.update_lambda(delta_lambda_A, delta_lambda_B); 
    
    % update iteration. 
    iter = iter + 1; 

    % update errors. 
    s_vector = beam.compute_s(s, C, gamma_ref, omega_ref); 
    f_int = beam.compute_f_int(n_gauss_points, C, gamma_ref, omega_ref);
    f_H1 = beam.compute_f_H1();
    f_H2 = beam.compute_f_H2();
    f_H = [f_H1; f_H2]; 
    residual = f_int + f_H - f_ext; 
    residual = residual(free_dof_res); 
    h_q = beam.compute_h_q();
    h_q = h_q(free_dof_h);
    error_h = MSE_error(h_q); 
    error_res = MSE_error(residual); 

end

%% plot beam. 
x_lim = 5; 
y_lim = 3; 
z_lim = 3;
scaling = 0;
title = "Force at node 2: [Fx, Fy, Fz] = [" + num2str(Fx) + "N, " + num2str(Fy) + "N, " + num2str(Fx) + "N ]"; 
beam.show_config(x_lim, y_lim, z_lim, scaling, title);


