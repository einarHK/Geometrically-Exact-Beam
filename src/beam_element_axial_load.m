clear; 
clc; 

% Axial elongation - 4 nonlinear equations solver. 
%% CONSTANTS. 
L = 1; % initial length. 
x1 = transpose([0, 0, 0]); % start coordinate. 
x2 = transpose([L, 0, 0]); % end coordinate. 
d1 = transpose([0, 0, 1]); 
d2 = transpose([0, 1, 0]); 
d3 = transpose([1, 0, 0]); 

kinematic_dof = 12;
n_constraints = 6; 
lambda_A = zeros(n_constraints, 1); 
lambda_B = zeros(n_constraints, 1);
fixed_dof = 1; 

gamma_ref = transpose([0, 0, 1]); 
omega_ref = transpose([0, 0, 0]); 

% initial e and s guess. 
e0 = [gamma_ref; omega_ref]; 
s0 = zeros(6,1); 

beam = Beam_elem(x1, x2, L, d3, d2, d1, d3, d2, d1, lambda_A, lambda_B, kinematic_dof, n_constraints, fixed_dof, gamma_ref, omega_ref);

beam.e0 = zeros(6,1); 
beam.s0 = zeros(6,1);
%% SIMULATION. 
% material property matrix. 
C = eye(6); 

load_steps = 10; 
fx = 1; 
fy = 0; 
fz = 0;

% dimensions for the strain and stress vectors. 
e_dim = 6; 
s_dim = 6;

% indices for the external force. 
fx_indx = kinematic_dof + 1; 
fy_indx = kinematic_dof + 2; 
fz_indx = kinematic_dof + 3;

% total degrees of freedom - kinematic, constraints, e, s. 
dof = (2 * kinematic_dof + 2 * n_constraints + e_dim + s_dim);
% external force vector. 
f_ext = zeros(2 * kinematic_dof, 1);

% tolerance and max iterations. 
TOL = 1e-6;
max_iter = 50; 

% fixed dof. 
q_fixed = 1:12;  
lambda_fixed = 37:42;

fixed_dof = [q_fixed, lambda_fixed]; % fix the first node and its constraints. 
all_dof = 1:dof; 
free_dof = setdiff(all_dof, fixed_dof); % free dofs.  

% number of gauss points. 
n_gauss_points = 1; 

%% SIMULATION.
total_iter = 0; 
for i=1:load_steps
   f_ext(fx_indx) = (i/load_steps) * fx;  
   f_ext(fy_indx) = (i/load_steps) * fy; 
   f_ext(fz_indx) = (i/load_steps) * fz; 
   
   iter = 0; 
   % Newtons method. 
   while (iter < max_iter)
        iter = iter + 1; % iteration counter. 
        delta_u = zeros(dof, 1); % solution vector. 
        
        % S matrix. 
        S = beam.compute_full_S_mat(0, n_gauss_points, C, gamma_ref, omega_ref);
        S = S(free_dof, free_dof); 
        
        % rhs vector. 
        f_int = beam.compute_f_int(n_gauss_points, C, gamma_ref, omega_ref); 
        f_H1 = beam.compute_f_H1(); 
        f_H2 = beam.compute_f_H2();
        f_H = [f_H1; f_H2]; 
        residual = f_int + f_H - f_ext; 
        hq = beam.compute_h_q(); 
        e_q = beam.compute_e(0, gamma_ref, omega_ref); 
        f = beam.e0 - e_q; % e - e(q)
        g = beam.s0 - C * beam.e0; % g(e,s)
        rhs = [-residual; -f;  -g; -hq]; % rhs vector. 
        rhs = rhs(free_dof); 

        % display norm and iterations. 
        disp("Iter: " + num2str(iter) + " Error: " + num2str(norm(rhs)));

        if (norm(rhs) < TOL)
            break; 
        end

        % solve for u. 
        delta_u(free_dof) = S \ rhs;
        dq = delta_u(1:2 * kinematic_dof); 
        de = delta_u(2 * kinematic_dof + 1: 2 * kinematic_dof + e_dim); 
        ds = delta_u(2 * kinematic_dof + e_dim + 1 : 2 * kinematic_dof + e_dim + s_dim); 
        dLambda = delta_u(2 * kinematic_dof + e_dim + s_dim + 1 : dof); 
        
        % update e. 
        beam.e0 = beam.e0 + de; 
        % update coordinates. 
        beam.update_params(dq(1:3), dq(13:15), dq(4:6), dq(7:9), dq(10:12), dq(16:18), dq(19:21), dq(22:24)); 
        % update s. 
        beam.s0 = beam.s0 + ds;     
        % update lambdas. 
        beam.update_lambda(dLambda(1:n_constraints), dLambda(n_constraints + 1: 2 * n_constraints)); 
   end
   total_iter = total_iter + iter; 
end

avg_iter = total_iter / load_steps; 
output = "Average number of iterations: " + avg_iter; 
disp(output); 
%% Display Strain and Stress. 
for i = 1:6
    stress = beam.s0(i); 
    strain = beam.e0(i); 
    output = "Strain component " + num2str(i) + " = " + num2str(strain); 
    output2 = "Stress component " + num2str(i) + " = " + num2str(stress); 
    disp(output); 
    disp(output2); 
    fprintf("\n"); 
end

% Calculate numerical displacement. 
dX = beam.x2_t - beam.x2; 
% dX_size = sqrt(dot(dX, dX)); 

% Calculate analytical displacement. 
dx_analytical = fx * L; % A=1, C=1 
dy_analytical = fy * L; % A=1, C=1
dz_analytical = fz * L; % A=1, C=1

output = "Numerical beam displacement : [Δx, Δy, Δz] = [ " + num2str(dX(1)) + "m, " + num2str(dX(2)) + "m, " + num2str(dX(3)) + "m ]";  
output2 = "Analytical beam displacement: [Δx, Δy, Δz] = [ " + num2str(dx_analytical) + "m, " + num2str(dy_analytical) + "m, " + num2str(dz_analytical) + "m ]";

disp(output); 
disp(output2);
%% PLOTS. 
% 
% e_min = -10; 
% e_max = 10; 
% C = 1; 
% N_points = 100; 
% noise = 0.1
% [e,s] = compute_linear_stress_strain(N_points, 1, e_min, e_max, noise); 

x_lim1 = -1; 
x_lim2 = fx * L + 2; 
y_lim1 = -1; 
y_lim2 = 1; 
z_lim1 = 0; 
z_lim2 = 2; 
scaling = 0; 
title_str = "Beam deformed config - F = [ " + num2str(fx) + "N, " + num2str(fy) + "N, " + num2str(fz) + "N ]"; 
scale_factor = 0.2; 

beam.show_config(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scaling, title_str, scale_factor)
