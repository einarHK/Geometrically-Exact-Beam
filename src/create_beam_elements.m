
% creates an array of beam elements, based on the coordinates input for
% each node. 
function beam_elems = create_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas)
    beam_elems = []; % initialize array. 
    % iterate over each beam element. 
    for i=1:n_elems
        x1 = beam_coordinates(:,i); % node 1.
        x2 = beam_coordinates(:,i+1); % node 2. 
        d1 = beam_directors(:, i); % directors at node 1. 
        d11 = d1(1:3); % shear director  - node 1.
        d12 = d1(4:6); 
        d13 = d1(7:9); % axis director  - node 1.
        d2 = beam_directors(:,i+1); % directors at node 2. 
        d21 = d2(1:3); % shear director - node 2. 
        d22 = d2(4:6); 
        d23 = d2(7:9); % axis director - node 2. 
        L0 = beam_lengths(i); % the initial length for the beam. 
        % assume unstressed reference configuration. 
        gamma_ref = [0; 0; 1]; 
        omega_ref = [0; 0; 0]; 
        % number of constraints on the element. 
        n_constraints = beam_constraints(i); 
        % free dof - kinematic variables. 
        dof = beam_dofs(i);
        % fixed dof for given element. 
        fixed_dof = beam_fixed_dofs(i); 
        % lambda vector - node 1. 
        lambda_A = lambdas(:,i); 
        % lambda vector - node 2. 
        lambda_B = lambdas(:,i+1);
        % create the beam class object. 
        beam = Beam_elem(x1, x2, L0, d13, d12, d11, d23, d22, d21, lambda_A, lambda_B, dof, n_constraints, fixed_dof, gamma_ref, omega_ref); 
        % add beam element to the array of beam elements. 
        beam_elems = [beam_elems, beam]; 
    end

end