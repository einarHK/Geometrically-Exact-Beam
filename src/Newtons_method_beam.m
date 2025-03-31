
% Newtons method - static case 
% Beam with multiple beam elements. 
function [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, TOL, f_ext, verbose, free_dof)
    % iteration counter.
    iter = 0;
    % main loop. 
    while (iter < max_iter)
        % increment iteration counter. 
        iter = iter + 1; 
        
        % initialize solution vector. 
        delta_u_lambda = zeros(beam.dof_per_node * beam.n_nodes + beam.constraint_per_node * beam.n_nodes, 1); 
        
        % compute S matrix. 
        S = beam.compute_S_mat(n_gauss_points, C); 
        S = S(free_dof, free_dof); 

        % compute the force terms. 
        % f_int internal force vector. 
        f_int = beam.compute_f_int(n_gauss_points, C); 
        % f_H force coming from the constraints. 
        f_H = beam.compute_f_H();
        
        % compute rhs vector. 
        % force residual vector. 
        residual = f_int + f_H - f_ext; 
        % constraint vector.    
        h_q = beam.compute_h_q();
        % create rhs vector. 
        rhs = [-residual; -h_q]; 
        rhs = rhs(free_dof); 

        % solve the linearized system of equations. 
        delta_u_lambda(free_dof) = S \ rhs; 

        % upate beam parameters. 
        beam.update_beam_params(delta_u_lambda); 

        % compute the norm of the rhs vector. 
        norm_rhs = norm(rhs); 

        % if verbose, show error at each step. 
        if (verbose ~= 0)
            output = "Error norm: " + num2str(norm_rhs) + ", Iterations: " + num2str(iter);
            disp([output]); 
        end

        % if norm of rhs less than tolerance, break out of loop. 
        if (norm_rhs < TOL)
            break;
        end

    end


end