
clear; 
clc; 
% Axial elongation - 4 nonlinear equations solver. 
%% CONSTANTS. 
n_elems = 1; % number of beam elements. 
L = 2; % initial length. 
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

% init e and s. 
e0 = zeros(6 * (n_elems + 1), 1); 
s0 = zeros(6 * (n_elems + 1), 1); 

beam.init_strain_stress(e0, s0);
%% SIMULATION. 
% material property matrix. 
C = eye(6); 

load_steps = 10; 
fx = 0; 
fy = 1; 
fz = 0;

% dimensions for the strain and stress vectors. 
e_dim = 6; 
s_dim = 6;

% indices for the external force. 
indx = n_elems + 1; 
fx_indx = (indx - 1) * beam.dof_per_node + 1; 
fy_indx = (indx - 1) * beam.dof_per_node + 2; 
fz_indx = (indx - 1) * beam.dof_per_node + 3;

% total degrees of freedom - kinematic, constraints, e, s. 
dof = (beam.n_nodes * beam.dof_per_node + beam.n_nodes * beam.constraint_per_node + beam.n_elements * e_dim + beam.n_elements * s_dim);
% external force vector. 
f_ext = zeros(beam.n_nodes * beam.dof_per_node, 1);

% tolerance and max iterations. 
TOL = 1e-6;
max_iter = 50; 

% fixed dof. 
fixed_nodes = [1];
q_fixed = []; 
lambda_fixed = []; 
fixed_dof = []; 
for j=1:length(fixed_nodes)
    node_indx = fixed_nodes(j); 
    q_vals = (node_indx - 1) * beam.dof_per_node + 1:node_indx * beam.dof_per_node;
    lambda_indx_start = beam.dof_per_node * beam.n_nodes + beam.n_elements * 6 + beam.n_elements * 6;  
    lambda_vals = (node_indx - 1) * beam.constraint_per_node + lambda_indx_start + 1: lambda_indx_start + node_indx * beam.constraint_per_node; 
    fixed_dof = [fixed_dof, q_vals]; 
    fixed_dof = [fixed_dof, lambda_vals];
end

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
        S = beam.compute_full_S_mat(n_gauss_points, C);
        S = S(free_dof, free_dof); 
        
        % rhs vector. 
        f_int = beam.compute_f_int(n_gauss_points, C); 
        f_H = beam.compute_f_H(); 
        residual = f_int + f_H - f_ext; 
        hq = beam.compute_h_q(); 
        e_q = beam.compute_strain(); 
        e = beam.get_strain_vector();
        f = e - e_q; % e - e(q)
        s = beam.get_stress_vector(); 
        g = beam.compute_g(e, s, C); % g(e,s) = s - C*e
        
        rhs = [-residual; -f;  -g; -hq]; % rhs vector. 
        rhs = rhs(free_dof); 
        
        % display norm and iterations. 
        disp("Iter: " + num2str(iter) + " Error: " + num2str(norm(rhs))); 
       

        if (norm(rhs) < TOL)
            break; 
        end

        % solve for u. 
        delta_u(free_dof) = S \ rhs;
        dq = delta_u(1:beam.n_nodes * beam.dof_per_node); 
        de = delta_u(beam.n_nodes * beam.dof_per_node + 1: beam.n_nodes * beam.dof_per_node + beam.n_elements * e_dim); 
        ds = delta_u(beam.n_nodes * beam.dof_per_node + beam.n_elements * e_dim + 1 : beam.n_nodes * beam.dof_per_node + beam.n_elements * (e_dim + s_dim)); 
        dLambda = delta_u(beam.n_nodes * beam.dof_per_node + beam.n_elements *( e_dim + s_dim) + 1 : dof); 
        

        % update e and s. 
        beam.update_beam_strain_stress(de, ds); 
        % update coordinates. 
        beam.update_beam_q(dq); 
        % update lambdas. 
        beam.update_beam_lambdas(dLambda); 
   end
   total_iter = total_iter + iter; 
end

avg_iter = total_iter / load_steps; 
output = "Average number of iterations: " + avg_iter; 
disp(output); 
%% Display Strain and Stress. 

% Calculate numerical displacement. 
beam_element = beam.beam_elements(beam.n_elements); 
dX = beam_element.x2_t - beam_element.x2; 
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
x_lim2 = fx + 3; 
y_lim1 = -1; 
y_lim2 = 5; 
z_lim1 = 0; 
z_lim2 = 3; 

scaling = 1; 
scale_factor = 0.5;
title = "Deformed beam, F = [" + num2str(fx) + "N, " + num2str(fy) + "N, " + num2str(fz) + "N]";

beam.plot_deformed_beam(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, 1, title);