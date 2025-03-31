% INPUT: number of beam elements n_elems, beam nodal coordinates 
% beam_coordinates, beam nodal directors beam_directors, beam initial
% lengths beam_lengths, beam constraints beam_constraints, beam kinematic
% degrees of freedom beam_dofs, beam fixed dof per node beam_fixed_dofs, beam node lagrange multipliers lambdas. 
% OUTPUT: array of beam elements of class type Beam_elem beam_elems. 
function beam_elems = create_curved_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas)
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
        % calculate gamma ref. 
        gamma_ref = compute_gamma(0, L0, x1, x2, d11, d12, d13, d21, d22, d23); 
        % calculate omega ref. 
        omega_ref = compute_omega(0, L0, d11, d12, d13, d21, d22, d23);
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