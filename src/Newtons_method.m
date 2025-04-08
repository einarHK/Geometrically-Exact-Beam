% Newtons method function - single beam element. 
% returns total number of iteration for convergence. 

function [iter] = Newtons_method(beam, n_gauss_points, C, gamma_ref, omega_ref, max_iter, TOL, free_dof, f_ext, n_constraints, verbose)
    iter = 0; 
    % main loop. 
    while (iter < max_iter)
        % increment iteration count. 
        iter = iter + 1; 
        
        % initialize solution vector. 
        delta_u_lambda = zeros(2 * beam.dof + 2 * beam.n_constraints, 1);
        
        % compute S matrix. 
        S = beam.compute_S_mat(n_gauss_points, C, gamma_ref, omega_ref); 
        S = S(free_dof, free_dof); 
        % disp(cond(S));
      
        % compute the force terms. 
        f_int = beam.compute_f_int(n_gauss_points, C, gamma_ref, omega_ref); % internal force term. 
        f_H1 = beam.compute_f_H1(); % force coming from node 1 constraints. 
        f_H2 = beam.compute_f_H2(); % forces coming from node 2 constraints. 
        f_H = [f_H1; f_H2]; 
        residual = f_int + f_H - f_ext; % residual force vector.

        h_q = beam.compute_h_q(); % constraint vector. 
        rhs = [-residual; -h_q]; 
        rhs = rhs(free_dof); % right hand side vector. 

        % solve linearized system of equations. 
        delta_u_lambda(free_dof) = S \ rhs; 
        
        % update beam parameters. 
        beam.update_config(delta_u_lambda);
 
        % extract the beam value increments. 
        % dx1 = delta_u_lambda(1:3); 
        % delta_d1_A = delta_u_lambda(4:6); 
        % delta_d2_A = delta_u_lambda(7:9); 
        % delta_d3_A = delta_u_lambda(10:12); 
        % dx2 = delta_u_lambda(13:15); 
        % delta_d1_B = delta_u_lambda(16:18); 
        % delta_d2_B = delta_u_lambda(19:21); 
        % delta_d3_B = delta_u_lambda(22:24); 
        % delta_lambda_A = delta_u_lambda(25:24 + n_constraints, 1); 
        % delta_lambda_B = delta_u_lambda(25 + n_constraints: 24 + 2*n_constraints); 
        
        % updating the beam parameters. 
        % beam.update_params(dx1, dx2, delta_d1_A, delta_d2_A, delta_d3_A, delta_d1_B, delta_d2_B, delta_d3_B);
        % beam.update_lambda(delta_lambda_A, delta_lambda_B); 
        
        % beam.display_node_pos();
        % compute the norm of the rhs vector. 
        norm_error = norm(rhs);
        
        % if verbose, show error at each step. 
        if (verbose ~= 0)
            output = "Error norm: " + num2str(norm_error) + ", Iterations: " + num2str(iter);
            % disp(delta_u_lambda)
            disp([output]); 
        end

        % if the norm less than tolerance, break out of loop. 
        if (norm_error < TOL)
            break;
        end
    end

end