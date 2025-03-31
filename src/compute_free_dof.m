
% returns the free dof indices, given array of fixed nodes. 
function [free_dof] = compute_free_dof(fixed_nodes, kinematic_dof, constraints_dof, n_nodes)
    all_dof = 1:kinematic_dof * n_nodes + constraints_dof * n_nodes;
    fixed_constraints_dof = []; 
    fixed_kinematic_dof = []; 
    for i=1:length(fixed_nodes)
        node = fixed_nodes(i); 
        j_start = (node - 1) * kinematic_dof + 1; 
        j_end = node * kinematic_dof; 
        k_start = kinematic_dof * n_nodes + (node - 1) * constraints_dof + 1; 
        k_end = kinematic_dof * n_nodes + node * constraints_dof; 
        for j=j_start:j_end
            fixed_kinematic_dof = [fixed_kinematic_dof, j];
        end
        for k=k_start:k_end
            fixed_constraints_dof = [fixed_constraints_dof, k];
        end
    end
    free_dof = setdiff(all_dof, [fixed_kinematic_dof, fixed_constraints_dof]);
end

