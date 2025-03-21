
% Multiple beam element structure. 
classdef Beam < handle
    properties
        n_elems; % total number of beam elements. 
        n_nodes; % total number of beam nodes. 
        beam_elements; % array of the beam elements. 
        dof_per_node; % number of degrees of freedom kinematic variables at each node.
        n_constraints; % number of constraints per node. 
    end

    methods
        % constructor. 
        function obj = Beam(n_elems, beam_elements, dof_per_node, n_constraints)
            obj.n_elems = n_elems;
            obj.n_nodes = n_elems + 1; 
            obj.beam_elements = beam_elements; 
            obj.dof_per_node = dof_per_node;
            obj.n_constraints = n_constraints; 
        end
        
        % computes the internal force vector for the beam, 
        % adding all the contributions from each beam element. 
        function f_int = compute_f_int(obj, n_gauss_points, C)
            % initialize the force vector. 
            f_int = zeros(obj.dof_per_node * obj.n_nodes, 1); 
            % iterate over each beam element. 
            for i=1:obj.n_elems
                i1 = 1 + (i - 1) * obj.dof_per_node ; 
                i2 = i1 + 2 * obj.dof_per_node - 1;
                beam = obj.beam_elements(i); 
                f_int_elem = beam.compute_f_int(n_gauss_points, C, beam.gamma_ref, beam.omega_ref);
                disp(f_int_elem)
                f_int(i1:i2) = f_int(i1:i2) + f_int_elem; % add to the total internal force. 
            end
        end
        
        % computes the force coming from the constraints. 
        function f_H = compute_f_H(obj)
            % initialize the force vector. 
            f_H = zeros(obj.dof_per_node * obj.n_nodes, 1); 
            for i=1:obj.n_elems
                i1 = (i-1) * obj.dof_per_node + 1; 
                i2 = (i) * obj.dof_per_node; 
                beam = obj.beam_elements(i); 
                f_H1 = beam.compute_f_H1(); % force coming from constraint at node 1. 
                f_H(i1:i2) = f_H(i1:i2) + f_H1; 
            end
            % add force coming from last constraint at end node. 
            f_H2 = beam.compute_f_H2(); 
            f_H(i2+1:obj.dof_per_node * obj.n_nodes) = f_H(i2+1:obj.dof_per_node * obj.n_nodes) + f_H2;
        end

        % plots the current configuration of the beam. 
        function show_config(obj, x_lim, y_lim, z_lim, scaling, title)
            % iterate over each beam element. 
            for i=1:obj.n_elems
                beam = obj.beam_elements(i);
                beam.show_config(x_lim, y_lim, z_lim, scaling, title)
            end
        end

    end

end




