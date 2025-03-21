% Newtons method function. 
% returns total number of iteration for convergence. 

function [iter] = Newtons_method(beam, s, C, gamma_ref, omega_ref, max_iter, TOL_q, TOL_lambda, f_ext, n_gauss_points)

% compute initial q. 
q = beam.compute_q(s); 
% initial stress & moments.
s_vector = beam.compute_N(s, C, gamma_ref, omega_ref); 
% iteration counter. 
iter = 0; 
% compute free_dof. 
free_dof = beam.compute_free_dof();

% compute initial error. 
h_q = beam.compute_h_nodes(); 
dof_h_q = length(h_q);

fy = Expand_vector(free_dof, dof_h_q); 
fx = Expand_vector(free_dof, dof_h_q);

error_h = MSE_error(h_q);
% compute residual. 
f_int = beam.compute_f_int(C, n_gauss_points, gamma_ref, omega_ref); 
H_h = beam.compute_H_h(); 
lambda = [beam.lambda_At; beam.lambda_Bt]; 
f_H = H_h' * lambda; 
residual = f_int + f_H - f_ext;
residual = residual(free_dof);
error_res = MSE_error(residual); 

while (((error_res > TOL_q) && (error_h > TOL_lambda)) || (iter < max_iter))
    % init solution vector. 
    delta_u = zeros(beam.dof * 2, 1); % u vector - displacement + directors. 
    delta_lambda = zeros(beam.n_constraints * 2, 1); % lambda vector - number of constraints for element. 
    delta_u_lambda = [delta_u; delta_lambda]; % combined vector - u + lambda.  
    
    % compute the S matrix - containing the tangent stiffness matrices. 
    S = beam.compute_S_matrix_h(s, C, s_vector);
    S = S(fy, fx); % index matrix S to avoid fixed indices. 
    
    % compute the residual vector.
    f_int = beam.compute_f_int(C, n_gauss_points, gamma_ref, omega_ref); 
    H_h = beam.compute_H_h(); 
    lambda = [beam.lambda_At; beam.lambda_Bt]; 
    f_H = H_h' * lambda; % force from constraints. 
    residual = f_int + f_H - f_ext;
    residual = residual(free_dof); 

    % compute the h(q) constraint vector. 
    h_q = beam.compute_h_nodes(); 

    % collect the residuals in the rhs vector. 
    rhs = [-residual; h_q];  

    % solve the linearized system of equations. 
    delta_u_lambda(fy) = S \ rhs; 
    i1 = 1; 
    i2 = 2 * beam.dof; 
    i3 = length(delta_u_lambda); 
    delta_u = delta_u_lambda(i1:i2); 
    delta_lambda = delta_u_lambda(i2+1:i3);
    % update beam values. 
    beam.update_params(delta_u); 
    beam.update_lagrange_coeffs(delta_lambda(1:beam.n_constraints), delta_lambda(beam.n_constraints+1:2*beam.n_constraints));
    
    % compute new stress vector. 
    s_vector = beam.compute_N(s, C, gamma_ref, omega_ref); 

    % update iteration counter. 
    iter = iter + 1; 

    % update errors. 
    f_int = beam.compute_f_int(C, n_gauss_points, gamma_ref, omega_ref);
    H_h = beam.compute_H_h(); 
    lambda = [beam.lambda_At; beam.lambda_Bt]; 
    f_H = H_h' * lambda;
    residual = f_int + f_H - f_ext; 
    residual = residual(free_dof);
    h_q = beam.compute_h_nodes(); 
    error_h = MSE_error(h_q); 
    error_res = MSE_error(residual); 
    
end

end

