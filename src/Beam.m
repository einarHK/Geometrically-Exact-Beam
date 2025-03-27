
% Beam class - consists of multiple beam elements. 
classdef Beam < handle
    properties
        beam_elements; % beam elements. 
        n_elements; % number of beam elements.  
        dof_per_node; % kinematic dof per node. 
        constraint_per_node; % constraints per node. 
        n_nodes; % total number of beam nodes. 
    end
    methods
        % constructor. 
        function obj = Beam(n_elements, beam_elements, dof_per_node, constraint_per_node)
            obj.beam_elements = beam_elements;
            obj.n_elements = n_elements;
            obj.dof_per_node = dof_per_node; 
            obj.constraint_per_node = constraint_per_node;
            obj.n_nodes = n_elements + 1;
        end
        
        % computes the free dof indices. 
        function free_dof = compute_free_dof(obj)
            fixed_kinematic_dof = []; 
            fixed_constraints = []; 
            % compute the fixed kinematic dof. 
            i2 = (obj.n_elements + 1) * obj.dof_per_node; % start index - constraints. 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                beam_fixed_dof = beam.fixed_dof; 
                if (beam_fixed_dof > 0)
                    i_start = 1 + (i - 1) * obj.dof_per_node; 
                    for j=i_start:beam_fixed_dof + i_start - 1
                        fixed_kinematic_dof = [fixed_kinematic_dof, j]; 
                    end
                    k_start = i2 + 1 + (i - 1) * obj.constraint_per_node; 
                    for k=k_start:k_start + obj.constraint_per_node - 1
                        fixed_constraints = [fixed_constraints, k];
                    end
                end
            end
            n_dof = (obj.n_elements + 1) * obj.dof_per_node + (obj.n_elements + 1) * obj.constraint_per_node;
            free_dof = 1:n_dof;
            free_dof = setdiff(free_dof, [fixed_kinematic_dof, fixed_constraints]); 
        end

        % internal force vector. 
        function f_int = compute_f_int(obj, n_gauss_points, C)
            f_int = zeros(obj.dof_per_node * (obj.n_elements + 1), 1); 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                f_contrib = beam.compute_f_int(n_gauss_points, C, beam.gamma_ref, beam.omega_ref); 
                i1 = 1 + (i - 1) * obj.dof_per_node;
                i2 = i1 - 1 + 2 * obj.dof_per_node; 
                f_int(i1:i2) = f_int(i1:i2) + f_contrib; 
            end
        end

        % force coming from constraints. 
        function f_H = compute_f_H(obj)
            f_H = zeros(obj.dof_per_node * (obj.n_elements + 1), 1);
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                f_H1 = beam.compute_f_H1();
                i1 = 1 + (i - 1) * obj.dof_per_node; 
                i2 = i1 + obj.dof_per_node - 1; 
                f_H(i1:i2) = f_H(i1:i2) + f_H1; 
            end
            i1 = 1 + obj.n_elements * obj.dof_per_node; 
            i2 = i1 + obj.dof_per_node - 1; 
            f_H2 = beam.compute_f_H2();
            f_H(i1:i2) = f_H(i1:i2) + f_H2; 
        end
        
        % computes the stiffness matrix - material and geometric.
        function Kt = compute_Kt(obj, n_gauss_points, C)
            % init Kt.
            Kt = zeros(obj.dof_per_node * (obj.n_elements + 1), obj.dof_per_node * (obj.n_elements + 1));
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                kt = beam.compute_Kt(n_gauss_points, C, beam.gamma_ref, beam.omega_ref);
                i1 = 1 + (i - 1) * beam.dof;
                i2 = i1 - 1 + 2 * beam.dof;
                Kt(i1:i2, i1:i2) = Kt(i1:i2, i1:i2) + kt; 
            end
        end
        
        % computes the dH_dq matrix. 
        function dH_dq = compute_dH_dq(obj)
            dH_dq = zeros(obj.dof_per_node * (obj.n_elements + 1), obj.dof_per_node * (obj.n_elements + 1));
            Z12 = zeros(obj.dof_per_node, obj.dof_per_node); 
            for i = 1:obj.n_elements
                beam = obj.beam_elements(i); 
                V1 = beam.compute_V1();
                dh_dq = [V1, Z12; Z12, Z12]; 
                i1 = 1 + (i - 1) * obj.dof_per_node; 
                i2 = i1 - 1 + 2 * obj.dof_per_node; 
                dH_dq(i1:i2, i1:i2) = dH_dq(i1:i2, i1:i2) + dh_dq;
            end
            V2 = beam.compute_V2(); % second node, last element. 
            dh_dq = [Z12, Z12; Z12, V2]; 
            i1 = obj.dof_per_node * (obj.n_elements - 1) + 1; 
            i2 = obj.dof_per_node * (obj.n_elements + 1); 
            dH_dq(i1:i2, i1:i2) =  dH_dq(i1:i2, i1:i2) + dh_dq; 

        end
        
        % computes the vector h=h(q).
        function h_q = compute_h_q(obj)
            h_q = zeros((obj.n_elements + 1) * obj.constraint_per_node, 1); 
            for i=1:obj.n_elements
                i1 = 1 + (i - 1) * obj.constraint_per_node; 
                i2 = i1 - 1 + obj.constraint_per_node;
                beam = obj.beam_elements(i);
                h1 = beam.compute_h1(); 
                h_q(i1:i2) = h_q(i1:i2) + h1;
            end
            h2 = beam.compute_h2();
            i2 = (1 + obj.n_elements) * obj.constraint_per_node; 
            i1 = i2 - obj.constraint_per_node + 1; 
            h_q(i1:i2) = h_q(i1:i2) + h2; 
        end

        % computes the derivative of h=h(q). 
        function H = compute_H(obj)
            H = zeros((obj.n_elements + 1) * obj.constraint_per_node, (obj.n_elements + 1) * obj.dof_per_node);
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                H1 = beam.compute_H_node1(); % H matrix node 1.  
                [n_rows, n_cols] = size(H1); % [6, 12] 
                Z = zeros(n_rows, n_cols); % zero matrix - dimensions 12x24. 
                H1_mat = [H1, Z; Z, Z];
                i1_1 = 1 + (i-1) * obj.constraint_per_node; 
                i1_2 = i1_1 - 1 + 2 * obj.constraint_per_node;
                i2_1 = 1 + (i-1) * obj.dof_per_node; 
                i2_2 = i2_1 - 1 + 2 * obj.dof_per_node;  
                H(i1_1:i1_2, i2_1:i2_2) = H(i1_1:i1_2, i2_1:i2_2) + H1_mat;
            end
            % second node, last beam. 
            i1_2 = (obj.n_elements + 1) * obj.constraint_per_node; 
            i1_1 = i1_2 - 2 * obj.constraint_per_node + 1; 
            i2_2 = (obj.n_elements + 1) * obj.dof_per_node; 
            i2_1 = i2_2 - 2 * obj.dof_per_node + 1; 
            H2 = beam.compute_H_node2();
            [n_rows, n_cols] = size(H2); % [6, 12] 
            Z = zeros(n_rows, n_cols); % zero matrix. 
            H2_mat = [Z, Z; Z, H2]; % the last matrix - node 2 last beam element. 
            H(i1_1:i1_2, i2_1:i2_2) = H(i1_1:i1_2, i2_1:i2_2) + H2_mat; 
        end
        
        % computes the S matrix for the linearized equations. 
        function S = compute_S_mat(obj, n_gauss_points, C)
            n_rows = obj.dof_per_node * (obj.n_elements + 1) + obj.constraint_per_node * (obj.n_elements + 1); 
            n_cols = obj.dof_per_node * (obj.n_elements + 1) + obj.constraint_per_node * (obj.n_elements + 1); 
            S = zeros(n_rows, n_cols); % init S matrix. 
            % compute submatrices. 
            S11 = obj.compute_Kt(n_gauss_points, C) + obj.compute_dH_dq(); 
            H = obj.compute_H(); 
            S12 = transpose(H); 
            S21 = H; 
            n_dof = obj.dof_per_node * (obj.n_elements + 1); 
            n_constraints = obj.constraint_per_node * (obj.n_elements + 1); 
            S(1:n_dof, 1:n_dof) = S(1:n_dof, 1:n_dof) + S11; 
            S(1:n_dof, n_dof+1:n_dof + n_constraints) = S(1:n_dof, n_dof+1:n_dof + n_constraints) + S12; 
            S(n_dof + 1:n_dof + n_constraints, 1:n_dof) = S(n_dof + 1:n_dof + n_constraints, 1:n_dof) + S21; 
        end

        % displays end node coordinates. 
        function display_end_node_pos(obj)
            beam = obj.beam_elements(obj.n_elements); 
            x = beam.x2_t(1); 
            y = beam.x2_t(2); 
            z = beam.x2_t(3); 
            output = "End node coordinates: (" + num2str(x) + "m, " + num2str(y) + "m, " + num2str(z) + "m )";
            disp(output);
        end

        % displays end node directors. 
        function display_end_node_directors(obj)
            beam = obj.beam_elements(obj.n_elements);
            % directors at the end node of the last beam. 
            d1_Bt = beam.d1_Bt; 
            d2_Bt = beam.d2_Bt; 
            d3_Bt = beam.d3_Bt; 
            output1 = "End node directors:";
            output2 = "[" + num2str(d1_Bt(1)) + ", " + num2str(d1_Bt(2)) + ", " + num2str(d1_Bt(3)) + "]";
            output3 = "[" + num2str(d2_Bt(1)) + ", " + num2str(d2_Bt(2)) + ", " + num2str(d2_Bt(3)) + "]";
            output4 = "[" + num2str(d3_Bt(1)) + ", " + num2str(d3_Bt(2)) + ", " + num2str(d3_Bt(3)) + "]";
            disp(output1); 
            disp(output2); 
            disp(output3);
            disp(output4); 
        end

        % update the beam parameters. 
        function update_beam_params(obj, delta_u_lambdas)
            i1 = 1; % start index - dof kinematic variables. 
            j1 = obj.dof_per_node * obj.n_nodes + 1; % start index - constraints.
            for i=1:obj.n_elements
                % index beam element. 
                beam = obj.beam_elements(i);
                i_start = i1 + (i - 1) * obj.dof_per_node;
                i_end = i_start - 1 + 2 * obj.dof_per_node;
                j_start = j1 + (i - 1) * obj.constraint_per_node; 
                j_end = j_start - 1 + 2 * obj.constraint_per_node; 
                % kinematic displacement values. 
                delta_u = delta_u_lambdas(i_start:i_end); 
                % lambda values.
                delta_lambdas = delta_u_lambdas(j_start:j_end); 
                % solution vector for the beam element. 
                delta_u_lambda = [delta_u;delta_lambdas]; 
                % update parameters given the solution vector. 
                beam.update_config(delta_u_lambda); 
            end
        end

        % plot the beam. 
        function plot(obj, x_lim, y_lim, z_lim, scaling, title)
            % plot each beam element. 
            for i=1:length(obj.beam_elements)-1
                beam = obj.beam_elements(i);
                beam.show_config(x_lim, y_lim, z_lim, scaling, title);
            end
            title_str = "Beam configuration - N elements = " + num2str(obj.n_elements);
            beam = obj.beam_elements(i+1);
            beam.show_config(x_lim, y_lim, z_lim, scaling, title_str);
        end

    end

end


