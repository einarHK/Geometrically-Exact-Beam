
clear; 
clc; 
%% CONSTANTS. 
L = 2; % initial length. 
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

C = eye(6); 

% check determinant of S. 
S = beam.compute_KKT_S(C, C, 1, zeros(6,1), zeros(6,1)); 
disp("Determinant of S: "); 
disp(det(S));
%% Compute e_tilde, s_tilde - Syntethic data. 
n_points = 1; 
e_tilde_vals = []; 
s_tilde_vals = []; 

for i=1:n_points
    strain = transpose([0, 0, 1 + (i-1), 0, 0, 0]);
    e_tilde_vals = [e_tilde_vals, strain]; 
    stress = C * strain; 
    s_tilde_vals = [s_tilde_vals, stress];  
end
%% Simulation. 
TOL = 1e-6; % tolerance value. 
load_steps = 1; % number of load steps. 

% indices for the external force. 
fx_indx = kinematic_dof + 1; 
fy_indx = kinematic_dof + 2; 
fz_indx = kinematic_dof + 3;
fx = 0; 
fy = 0; 
fz = 0;

e_tilde_dim = 6; 
s_tilde_dim = 6; 
e_dim = 6; 
s_dim = 6; 
chi_dim = 24; 
mu_dim = 6; 
gamma_dim = 12; 
xi_dim = 6;

% total degrees of freedom - kinematic, constraints, e, s. 
dof = (e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + 2 * n_constraints + e_dim + s_dim + chi_dim + mu_dim + gamma_dim + xi_dim);
% external force vector. 
f_ext = zeros(2 * kinematic_dof, 1);

% fixed dof. 
q_fixed = e_tilde_dim + s_tilde_dim + 1: e_tilde_dim + s_tilde_dim + 12;  
lambda_fixed = e_tilde_dim + s_tilde_dim + 25:e_tilde_dim + s_tilde_dim + 30;

fixed_dof = [q_fixed, lambda_fixed]; % fix the first node and its constraints. 
all_dof = 1:dof; 
free_dof = setdiff(all_dof, fixed_dof); % free dofs.  

% number of gauss points. 
n_gauss_points = 1; 
max_iter = 50;  

%% 

% loop over each synthetic strain-stress value pair. 
for i=1:n_points
    % get synthetic strain-stress value pair. 
    e_tilde = e_tilde_vals(:,i); 
    s_tilde = s_tilde_vals(:,i);
    
    for j=1:load_steps
        f_ext(fx_indx) = (j/load_steps) * fx;  
        f_ext(fy_indx) = (j/load_steps) * fy; 
        f_ext(fz_indx) = (j/load_steps) * fz;

        iter = 0; 
        
        % Newtons method. 
        while (iter < max_iter)
            % init solution vector. 
            delta_x = zeros(dof, 1); 
            % Compute the variation of the Lagrangian. 
            dL = -beam.compute_dL(C, C, e_tilde, s_tilde, n_gauss_points, f_ext);  
            % S matrix. 
            S = beam.compute_KKT_S(C, C, n_gauss_points, e_tilde, s_tilde); 
            % solution vector. 
            delta_x(free_dof) = S(free_dof, free_dof) \ dL(free_dof); 

            if (norm(delta_x) < TOL)                
                break;
            end

            % update beam parameters. 
            delta_e_tilde = delta_x(1:e_tilde_dim); 
            delta_s_tilde = delta_x(e_tilde_dim + 1: e_tilde_dim + s_tilde_dim); 
            delta_q = delta_x(e_tilde_dim + s_tilde_dim + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof); 
            delta_e = delta_x(e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim);
            delta_s = delta_x(e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim);
            delta_lambda = delta_x(e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints); 
            delta_chi = delta_x(e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + chi_dim); 
            delta_mu = delta_x(e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + chi_dim + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + chi_dim + mu_dim);
            delta_gamma = delta_x(e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + chi_dim + mu_dim + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + chi_dim + mu_dim + gamma_dim);
            delta_xi = delta_x(e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + chi_dim + mu_dim + gamma_dim + 1: e_tilde_dim + s_tilde_dim + 2 * kinematic_dof + e_dim + s_dim + 2 * n_constraints + chi_dim + mu_dim + gamma_dim + xi_dim);
                        
            % update position and directors. 
            beam.update_params(delta_q(1:3), delta_q(13:15), delta_q(4:6), delta_q(7:9), delta_q(10:12), delta_q(16:18), delta_q(19:21), delta_q(22:24));
            % update lambda coeffs. 
            beam.update_lambda(delta_lambda(1:6), delta_lambda(7:12)); 
            % update e strain. 
            beam.e0 = beam.e0 + delta_e; 
            % update s stress. 
            beam.s0 = beam.s0 + delta_s; 
            % update chi coeffs. 
            beam.update_chi(delta_chi);
            % update mu coeffs. 
            beam.update_mu(delta_mu); 
            % update gamma coeffs. 
            beam.update_gamma(delta_gamma); 
            % update xi coeffs. 
            beam.update_xi(delta_xi); 
            
            iter = iter + 1; 

        end
    end

end



