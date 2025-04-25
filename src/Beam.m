
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

        % for initialization of e and s vectors. 
        function init_strain_stress(obj, e0, s0)            
            for i=1:obj.n_elements
                i_start = 1 + (i-1) * 6; 
                i_end = i_start + 5; 
                beam = obj.beam_elements(i); 
                beam.e0 = e0(i_start:i_end); 
                beam.s0 = s0(i_start:i_end); 
            end
        end
        
        % returns the beam coordinates in the undeformed state. 
        function [coordinates] = get_undeformed_beam_coords(obj)
            coordinates = []; 
            % iterate over each beam element.
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                x1 = beam.x1(1); 
                y1 = beam.x1(2); 
                z1 = beam.x1(3); 
                coordinate = transpose([x1, y1, z1]);
                coordinates = [coordinates, coordinate]; 
            end
            x2 = beam.x2(1); 
            y2 = beam.x2(2); 
            z2 = beam.x2(3); 
            coordinate = transpose([x2, y2, z2]);
            coordinates = [coordinates, coordinate]; 
        end

        % returns the deformed beam coordinates. 
        function [coordinates] = get_deformed_beam_coords(obj)
            coordinates = []; 
            % iterate over each beam element.
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                x1 = beam.x1_t(1); 
                y1 = beam.x1_t(2); 
                z1 = beam.x1_t(3); 
                coordinate = transpose([x1, y1, z1]);
                coordinates = [coordinates, coordinate]; 
            end
            x2 = beam.x2_t(1); 
            y2 = beam.x2_t(2); 
            z2 = beam.x2_t(3); 
            coordinate = transpose([x2, y2, z2]);
            coordinates = [coordinates, coordinate]; 
        end
        
        % returns the beam directors in the undeformed state. 
        function [directors] = get_undeformed_beam_directors(obj)
            directors = []; 
            % iterate over each beam element.
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                d1 = beam.d1_A; 
                d2 = beam.d2_A; 
                d3 = beam.d3_A; 
                d_vector = [d1;d2;d3]; 
                directors = [directors, d_vector]; 
            end
            d1 = beam.d1_B; 
            d2 = beam.d2_B; 
            d3 = beam.d3_B; 
            d_vector = [d1;d2;d3]; 
            directors = [directors, d_vector]; 
        end

        % returns the beam directors in the undeformed state. 
        function [directors] = get_deformed_beam_directors(obj)
            directors = []; 
            % iterate over each beam element.
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                d1 = beam.d1_At; 
                d2 = beam.d2_At; 
                d3 = beam.d3_At; 
                d_vector = [d1;d2;d3]; 
                directors = [directors, d_vector]; 
            end
            d1 = beam.d1_Bt; 
            d2 = beam.d2_Bt; 
            d3 = beam.d3_Bt; 
            d_vector = [d1;d2;d3]; 
            directors = [directors, d_vector]; 
        end

        % computes the kinematic fixed dof. 
        function fixed_kinematic_dof = compute_kinematic_dof(obj)
            fixed_kinematic_dof = [];
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                beam_fixed_dof = beam.fixed_dof; 
                if (beam_fixed_dof > 0)
                    j_start = (i - 1) * obj.dof_per_node;
                    for j=1:beam_fixed_dof
                        fixed_kinematic_dof = [fixed_kinematic_dof, j + j_start];
                    end
                end
            end
        end

        % computes the fixed constraints. 
        function fixed_constraints = compute_fixed_constraints(obj)
            fixed_constraints = []; 
            indx_start = obj.dof_per_node * obj.n_nodes; 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                if (beam.fixed_dof > 0)
                    for j=1:obj.constraint_per_node
                        j_indx = indx_start + i * obj.constraint_per_node + j; 
                        fixed_constraints = [fixed_constraints, j_indx];
                    end
                end
            end
        end

        % computes the free dof indices. 
        function free_dof = compute_free_dof(obj, fix_end_node)
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
          
            % free_dof = setdiff(free_dof, [fixed_kinematic_dof, fixed_constraints]); 
            if (fix_end_node > 0)
                end_node_dof = (obj.n_elements) * obj.dof_per_node+1:obj.n_nodes * obj.dof_per_node;
                end_node_constraints = (obj.n_nodes * obj.dof_per_node + obj.n_elements * obj.constraint_per_node + 1):n_dof;
                fixed_node_array = [fixed_kinematic_dof, fixed_constraints, end_node_dof, end_node_constraints];
                free_dof = setdiff(free_dof, fixed_node_array); 
            else
                free_dof = setdiff(free_dof, [fixed_kinematic_dof, fixed_constraints]); 
            end

        end

        % computes the strains for each beam, based on current
        % configuration - assume s = 0 for interpolating on the beam. 
        function e = compute_strain(obj)
            e = zeros(6 * obj.n_elements, 1); 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                e_q = beam.compute_e(0, beam.gamma_ref, beam.omega_ref); 
                i_start = 1 + (i - 1) * 6; 
                i_end = i_start + 5; 
                e(i_start:i_end) = e(i_start:i_end) + e_q;
            end
        end

        % computes the stresses for each beam, based on current 
        % configuration. 
        function s = compute_stress(obj, C)
            s = zeros(6 * obj.n_elements, 1);
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                s_val = beam.compute_s(0, C, beam.gamma_ref, beam.omega_ref); 
                i_start = 1 + (i - 1) * 6; 
                i_end = i_start + 5; 
                s(i_start:i_end) = s(i_start:i_end) + s_val;
            end
        end

        % returns current strain vector. 
        function e = get_strain_vector(obj)
            e = zeros(6 * obj.n_elements, 1); 
            for i = 1:obj.n_elements
                i_start = 1 + (i - 1) * 6; 
                i_end = i_start + 5; 
                beam = obj.beam_elements(i);
                e(i_start:i_end) =  e(i_start:i_end) + beam.e0;
            end
        end
        
        % returns current stress vector. 
        function s = get_stress_vector(obj)
            s = zeros(6 * obj.n_elements, 1); 
            for i = 1:obj.n_elements
                i_start = 1 + (i - 1) * 6; 
                i_end = i_start + 5; 
                beam = obj.beam_elements(i);
                s(i_start:i_end) =  s(i_start:i_end) + beam.s0;
            end
        end
        
        % computes g(e,s) = s - C * e. 
        function g = compute_g(obj, e, s, C)            
            g = zeros(6 * obj.n_elements, 1);
            for i=1:obj.n_elements
                i_start = 1 + (i - 1) * 6; 
                i_end = i_start + 5; 
                % beam = obj.beam_elements(i); 
                g(i_start:i_end) = s(i_start:i_end) - C * e(i_start:i_end); 
            end
        end

        % returns the stress and strain pair for the given elements. 
        function [stresses, strains] = compute_stress_strain(obj, elements, C)
            stresses = []; 
            strains = []; 
            for i=1:length(elements)
                element = elements(i); 
                beam_element = obj.beam_elements(element); 
                strain = beam_element.compute_e(0, beam_element.gamma_ref, beam_element.omega_ref); 
                stress = beam_element.compute_s(0, C, beam_element.gamma_ref, beam_element.omega_ref); 
                stresses = [stresses, stress]; 
                strains = [strains, strain]; 
            end
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
                % get the stiffness matrix from the beam element. 
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

        % The KKT matrix S_fix components. 
        % s - interpolation value. 
        function S_qq = compute_S_qq(obj, s)
            
            
        end

        function S_qs = compute_S_qs(obj, s)


        end


        % computes the full S matrix, when considering 4 nonlinear
        % equations instead of 2. 
        function S11 = compute_S11_component(obj, n_gauss_points, C)
            % init empty matrix. 
            S11 = zeros(obj.n_nodes * obj.dof_per_node, obj.n_nodes * obj.dof_per_node);
            for i = 1:obj.n_elements
                beam = obj.beam_elements(i); 
                S11_component = beam.compute_U2_integral_2(n_gauss_points, C, beam.gamma_ref, beam.omega_ref) + beam.compute_V_tot();
                i_start = obj.dof_per_node * (i - 1) + 1; 
                i_end = i_start - 1 + 2 * obj.dof_per_node; 
                S11(i_start:i_end, i_start:i_end) = S11(i_start:i_end, i_start:i_end) + S11_component;
            end            
        end
        
        function S12 = compute_S12_component(obj)
            S12 = zeros(obj.n_nodes * obj.dof_per_node, obj.n_elements * 6); 
        end

        function S13 = compute_S13_component(obj, n_gauss_points)
            S13 = zeros(obj.n_nodes * obj.dof_per_node, obj.n_elements * 6); 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                [weights, points] = beam.gauss_quadrature(n_gauss_points); 
                dx = beam.L0 / 2; 
                B_mat = zeros(2 * obj.dof_per_node, 6); 
                for j=1:n_gauss_points
                    point = points(j); 
                    weight = weights(j); 
                    B = beam.compute_B(point); 
                    B_mat = B_mat + transpose(B) * (weight * dx);  
                end
                i1 = 1 + (i - 1) * obj.dof_per_node; 
                i2 = i1 - 1 + 2 * obj.dof_per_node; 
                j1 = 1 + (i - 1) * 6; 
                j2 = j1 + 5;
                S13(i1:i2,j1:j2) = S13(i1:i2,j1:j2) + B_mat; 
            end
        end

        function S14 = compute_S14_component(obj)
            S14 = zeros(obj.n_nodes * obj.dof_per_node, obj.n_nodes * obj.constraint_per_node);     
            for i = 1:obj.n_elements
                beam = obj.beam_elements(i);
                H1 = beam.compute_H_node1(); 
                Z12 = zeros(obj.dof_per_node, obj.constraint_per_node);
                S14_component = [transpose(H1), Z12; Z12, Z12]; 
                i1 = 1 + (i - 1) * obj.dof_per_node; 
                i2 = i1 - 1 + 2 * obj.dof_per_node; 
                j1 = 1 + (i - 1) * obj.constraint_per_node; 
                j2 = j1 - 1 + 2 * obj.constraint_per_node;
                S14(i1:i2, j1:j2) = S14(i1:i2, j1:j2) + S14_component; 
            end
            H2 = beam.compute_H_node2();
            S14_component = [Z12, Z12; Z12, transpose(H2)]; 
            S14(i1:i2, j1:j2) = S14(i1:i2, j1:j2) + S14_component; 
        end

        function S21 = compute_S21_component(obj)
            S21 = zeros(6 * obj.n_elements, obj.n_nodes * obj.dof_per_node); 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                B = beam.compute_B(0);
                i1 = 1 + (i - 1) * 6; 
                i2 = i1 + 5; 
                j1 = 1 + obj.dof_per_node * (i - 1); 
                j2 = j1 - 1 + 2 * obj.dof_per_node; 
                S21(i1:i2, j1:j2) = S21(i1:i2, j1:j2) - B; 
            end
        end

        function S22 = compute_S22_component(obj)
            S22 = eye(6 * obj.n_elements); 
        end

        function S23 = compute_S23_component(obj) 
            S23 = zeros(6 * obj.n_elements, 6 * obj.n_elements); 
        end
        
        function S24 = compute_S24_component(obj)
           S24 = zeros(6 * obj.n_elements, obj.n_nodes * obj.constraint_per_node); 
        end

        function S31 = compute_S31_component(obj)
            S31 = zeros(6 * obj.n_elements, obj.n_nodes * obj.dof_per_node);
        end

        function S32 = compute_S32_component(obj, C)
            S32 = zeros(6 * obj.n_elements, 6 * obj.n_elements); 
            for i=1:obj.n_elements
                i1 = 1 + (i-1) * 6; 
                i2 = i1 + 5; 
                S32(i1:i2, i1:i2) = S32(i1:i2, i1:i2) -C; 
            end
        end
        
        function S33 = compute_S33_component(obj)
            S33 = eye(6 * obj.n_elements);
        end

        function S34 = compute_S34_component(obj)
            S34 = zeros(6 * obj.n_elements, obj.n_nodes * obj.constraint_per_node); 
        end

        function S41 = compute_S41_component(obj)
            S41 = zeros(obj.n_nodes * obj.constraint_per_node, obj.n_nodes * obj.dof_per_node); 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                H1 = beam.compute_H_node1(); 
                Z12 = zeros(obj.constraint_per_node, obj.dof_per_node); 
                S41_component = [H1, Z12; Z12, Z12]; 
                i1 = 1 + (i - 1) * obj.constraint_per_node; 
                i2 = i1 - 1 + 2 * obj.constraint_per_node; 
                j1 = 1 + (i - 1) * obj.dof_per_node; 
                j2 = j1 - 1 + 2 * obj.dof_per_node; 
                S41(i1:i2, j1:j2) = S41(i1:i2, j1:j2) + S41_component; 
            end
            H2 = beam.compute_H_node2(); 
            S41_component = [Z12, Z12; Z12, H2]; 
            S41(i1:i2, j1:j2) = S41(i1:i2, j1:j2) + S41_component; 
        end

        function S42 = compute_S42_component(obj)
            S42 = zeros(obj.n_nodes * obj.constraint_per_node, 6 * obj.n_elements);
        end

        function S43 = compute_S43_component(obj)
            S43 = zeros(obj.n_nodes * obj.constraint_per_node, 6 * obj.n_elements);
        end

        function S44 = compute_S44_component(obj)
            S44 = zeros(obj.n_nodes * obj.constraint_per_node, obj.n_nodes * obj.constraint_per_node); 
        end

        function S = compute_full_S_mat(obj, n_gauss_points, C)
            S11 = obj.compute_S11_component(n_gauss_points, C); 
            S12 = obj.compute_S12_component(); 
            S13 = obj.compute_S13_component(n_gauss_points);
            S14 = obj.compute_S14_component();
            S21 = obj.compute_S21_component(); 
            S22 = obj.compute_S22_component(); 
            S23 = obj.compute_S23_component(); 
            S24 = obj.compute_S24_component(); 
            S31 = obj.compute_S31_component(); 
            S32 = obj.compute_S32_component(C); 
            S33 = obj.compute_S33_component(); 
            S34 = obj.compute_S34_component(); 
            S41 = obj.compute_S41_component(); 
            S42 = obj.compute_S42_component(); 
            S43 = obj.compute_S43_component(); 
            S44 = obj.compute_S44_component();
            S = [S11, S12, S13, S14; S21, S22, S23, S24; S31, S32, S33, S34; S41, S42, S43, S44]; 
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

        % update the beam strain and stresses for each element. 
        function update_beam_strain_stress(obj, delta_e, delta_s)            
            for i=1:obj.n_elements
                i_start = 1 + 6 * (i - 1); 
                i_end = i_start + 5; 
                beam = obj.beam_elements(i);
                beam.e0 = beam.e0 + delta_e(i_start:i_end); 
                beam.s0 = beam.s0 + delta_s(i_start:i_end); 
            end
        end

        % update the beam positions. 
        function update_beam_q(obj, delta_q)           
            for i=1:obj.n_elements
                i_start = 1 + (i - 1) * obj.dof_per_node;
                i_end = i_start - 1 + 2 * obj.dof_per_node; 
                dq = delta_q(i_start:i_end); 
                dx1 = dq(1:3); 
                dx2 = dq(13:15); 
                delta_d1_A = dq(4:6); 
                delta_d2_A = dq(7:9); 
                delta_d3_A = dq(10:12); 
                delta_d1_B = dq(16:18); 
                delta_d2_B = dq(19:21); 
                delta_d3_B = dq(22:24); 
                beam = obj.beam_elements(i); 
                beam.update_params(dx1, dx2, delta_d1_A, delta_d2_A, delta_d3_A, delta_d1_B, delta_d2_B, delta_d3_B)
            end
        end
        
        % update beam lambda multipliers at each node. 
        function update_beam_lambdas(obj, delta_lambdas)
            for i=1:obj.n_elements
                i_start = 1 + (i - 1) * obj.constraint_per_node; 
                i1 = i_start - 1 + 1 * obj.constraint_per_node; 
                i2 = i_start - 1 + 2 * obj.constraint_per_node;
                beam = obj.beam_elements(i); 
                delta_lambda_A = delta_lambdas(i_start:i1); 
                delta_lambda_B = delta_lambdas(i1 + 1:i2); 
                beam.update_lambda(delta_lambda_A, delta_lambda_B); 
            end
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

        % returns the current coordinate values for each beam
        % element. 
        function [X, Y, Z] = get_beam_coords(obj)
            X = []; 
            Y = []; 
            Z = []; 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i); 
                x = beam.x1_t(1);
                y = beam.x1_t(2); 
                z = beam.x1_t(3);
                X = [X, x]; 
                Y = [Y, y]; 
                Z = [Z, z]; 
            end
            x2 = beam.x2_t(1);
            y2 = beam.x2_t(2); 
            z2 = beam.x2_t(3); 
            X = [X, x2]; 
            Y = [Y, y2]; 
            Z = [Z, z2];
        end
        
        % returns the current beam director values for each beam element. 
        function [directors] = get_beam_directors(obj)
            directors = []; 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                d1 = beam.d1_At; 
                d2 = beam.d2_At; 
                d3 = beam.d3_At; 
                d_vector = [d1;d2;d3]; 
                directors = [directors, d_vector];
            end 
        end
        
        % returns array of all nodal coordinates. 
        function [nodal_coordinates] = compute_nodal_coords(obj)
            nodal_coordinates = []; 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                x1 = beam.x1_t(1); 
                y1 = beam.x1_t(2); 
                z1 = beam.x1_t(3); 
                coord = transpose([x1, y1, z1]); 
                nodal_coordinates = [nodal_coordinates, coord]; 
            end
            x2 = beam.x2_t(1); 
            y2 = beam.x2_t(2); 
            z2 = beam.x2_t(3); 
            coord = transpose([x2, y2, z2]); 
            nodal_coordinates = [nodal_coordinates, coord]; 
        end

        % returns array of all the nodal directors. 
        function [directors] = compute_nodal_directors(obj)
            directors = []; 
            % iterate over each beam element. 
            for i=1:obj.n_elements
                beam = obj.beam_elements(i);
                d1 = beam.d1_At; 
                d2 = beam.d2_At; 
                d3 = beam.d3_At; 
                d_vector = [d1; d2; d3]; 
                directors = [directors, d_vector];
            end
            d1 = beam.d1_Bt; 
            d2 = beam.d2_Bt; 
            d3 = beam.d3_Bt; 
            d_vector = [d1; d2; d3]; 
            directors = [directors, d_vector]; 
        end
        
        % displays the current coordinate positions for the nodes given as
        % input. 
        function display_node_pos(obj, node_indices)
            % compute the list of all the nodal positions. 
            all_coords = obj.compute_nodal_coords();
            % iterate over each node index. 
            for i=1:length(node_indices)
                % get the node number.
                node_indx = node_indices(i);
                % get the coordinate. 
                coordinate = all_coords(:,node_indx);
                x1 = coordinate(1); 
                y1 = coordinate(2); 
                z1 = coordinate(3); 
                output = "Node " + num2str(node_indx) + " position: ( " + num2str(x1) + "m, " + num2str(y1) + ", "  + num2str(z1) + " )"; 
                disp(output);             
            end
        end

        % displays node directors. 
        function display_node_directors(obj, node_indices)
            all_directors = obj.compute_nodal_directors();
            for i=1:length(node_indices)
                node_indx = node_indices(i);
                d_vector = all_directors(:,node_indx); 
                d1 = d_vector(1:3); 
                d2 = d_vector(4:6);
                d3 = d_vector(7:9);                
                output = "Node " + num2str(node_indx) + " directors:"; 
                output2 = "d1 = (" + num2str(d1(1)) + ", " + num2str(d1(2)) + ", " + num2str(d1(3)) + " )";
                output3 = "d2 = (" + num2str(d2(1)) + ", " + num2str(d2(2)) + ", " + num2str(d2(3)) + " )"; 
                output4 = "d3 = (" + num2str(d3(1)) + ", " + num2str(d3(2)) + ", " + num2str(d3(3)) + " )"; 
                disp(output); 
                disp(output2);
                disp(output3); 
                disp(output4);
            end
        end
       
        % plot the beam. 
        function plot(obj, x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scaling, title, scale_factor)
            % plot each beam element. 
            for i=1:length(obj.beam_elements)-1
                beam = obj.beam_elements(i);
                beam.show_config(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scaling, title, scale_factor);
            end
            title_str = "Beam configuration - N elements = " + num2str(obj.n_elements);
            beam = obj.beam_elements(i+1);
            beam.show_config(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scaling, title_str, scale_factor);
        end
        
        % plots the undeformed state of the beam. 
        function plot_undeformed_state(obj, x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, title_str)
            beam_coordinates = obj.get_undeformed_beam_coords(); 
            beam_directors = obj.get_undeformed_beam_directors();
            for i=1:obj.n_elements
                p1 = beam_coordinates(:,i);
                p2 = beam_coordinates(:,i+1);
                x1 = p1(1); 
                y1 = p1(2); 
                z1 = p1(3); 
                X = [p1(1), p2(1)]; 
                Y = [p1(2), p2(2)]; 
                Z = [p1(3), p2(3)];
                d_vector = beam_directors(:,i); 
                d1 = d_vector(1:3); 
                d2 = d_vector(4:6); 
                d3 = d_vector(7:9);
                plot3(X, Y, Z, color="black", Marker="o", MarkerFaceColor="blue", MarkerSize=2);  
                hold on; 
                % plot directors. 
                d11_x = d1(1); 
                d11_y = d1(2); 
                d11_z = d1(3); 
                d12_x = d2(1); 
                d12_y = d2(2); 
                d12_z = d2(3); 
                d13_x = d3(1); 
                d13_y = d3(2); 
                d13_z = d3(3); 
                plot3([x1, x1 + scale * d11_x], [y1, y1 + scale *d11_y], [z1, z1 + scale*d11_z], Marker=">", LineStyle='-', MarkerSize=2, Color="blue");
                hold on;
                plot3([x1, x1 + scale*d12_x], [y1, y1 + scale*d12_y], [z1, z1 + scale*d12_z], Marker=">", LineStyle='-', MarkerSize=2, Color="blue");
                hold on;
                plot3([x1, x1 + scale*d13_x], [y1, y1 + scale*d13_y], [z1, z1 + scale*d13_z], Marker=">", LineStyle='-', MarkerSize=2, Color="blue");
                hold on;
            end
            x1 = p2(1); 
            y1 = p2(2); 
            z1 = p2(3); 
            d_vector = beam_directors(:,i+1); 
            d1 = d_vector(1:3); 
            d2 = d_vector(4:6); 
            d3 = d_vector(7:9);
            d11_x = d1(1); 
            d11_y = d1(2); 
            d11_z = d1(3); 
            d12_x = d2(1); 
            d12_y = d2(2); 
            d12_z = d2(3); 
            d13_x = d3(1); 
            d13_y = d3(2); 
            d13_z = d3(3); 
            plot3([x1, x1 + scale * d11_x], [y1, y1 + scale *d11_y], [z1, z1 + scale*d11_z], Marker=">", LineStyle='-', MarkerSize=2, Color="blue");
            hold on;
            plot3([x1, x1 + scale*d12_x], [y1, y1 + scale*d12_y], [z1, z1 + scale*d12_z], Marker=">", LineStyle='-', MarkerSize=2, Color="blue");
            hold on;
            plot3([x1, x1 + scale*d13_x], [y1, y1 + scale*d13_y], [z1, z1 + scale*d13_z], Marker=">", LineStyle='-', MarkerSize=2,  Color="blue");
            hold on;
            % set axis lim. 
            grid on;
            xlim([x_lim1, x_lim2]); 
            ylim([y_lim1, y_lim2]); 
            zlim([z_lim1, z_lim2]);
            xlabel("x-axis"); 
            ylabel("y-axis"); 
            zlabel("z-axis")
            % title of plot. 
            title([title_str]);
            % hold on; 
        end

        % plots the deformed beam coordinates. 
        function plot_deformed_beam(obj,x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, title_str)
            beam_coordinates = obj.get_deformed_beam_coords(); 
            beam_directors = obj.get_deformed_beam_directors();
            for i=1:obj.n_elements
                p1 = beam_coordinates(:,i);
                p2 = beam_coordinates(:,i+1);
                x1 = p1(1); 
                y1 = p1(2); 
                z1 = p1(3); 
                X = [p1(1), p2(1)]; 
                Y = [p1(2), p2(2)]; 
                Z = [p1(3), p2(3)];
                d_vector = beam_directors(:,i); 
                d1 = d_vector(1:3); 
                d2 = d_vector(4:6); 
                d3 = d_vector(7:9);
                plot3(X, Y, Z, color="red", Marker="o", MarkerFaceColor="blue", MarkerSize=2);  
                hold on; 
                % plot directors. 
                d11_x = d1(1); 
                d11_y = d1(2); 
                d11_z = d1(3); 
                d12_x = d2(1); 
                d12_y = d2(2); 
                d12_z = d2(3); 
                d13_x = d3(1); 
                d13_y = d3(2); 
                d13_z = d3(3); 
                plot3([x1, x1 + scale * d11_x], [y1, y1 + scale *d11_y], [z1, z1 + scale*d11_z], Marker=">", LineStyle='-', MarkerSize=2, Color="red");
                hold on;
                plot3([x1, x1 + scale*d12_x], [y1, y1 + scale*d12_y], [z1, z1 + scale*d12_z], Marker=">", LineStyle='-', MarkerSize=2, Color="red");
                hold on;
                plot3([x1, x1 + scale*d13_x], [y1, y1 + scale*d13_y], [z1, z1 + scale*d13_z], Marker=">", LineStyle='-', MarkerSize=2, Color="red");
                hold on;
            end
            x1 = p2(1); 
            y1 = p2(2); 
            z1 = p2(3); 
            d_vector = beam_directors(:,i+1); 
            d1 = d_vector(1:3); 
            d2 = d_vector(4:6); 
            d3 = d_vector(7:9);
            d11_x = d1(1); 
            d11_y = d1(2); 
            d11_z = d1(3); 
            d12_x = d2(1); 
            d12_y = d2(2); 
            d12_z = d2(3); 
            d13_x = d3(1); 
            d13_y = d3(2); 
            d13_z = d3(3); 
            plot3([x1, x1 + scale * d11_x], [y1, y1 + scale *d11_y], [z1, z1 + scale*d11_z], Marker=">", LineStyle='-', MarkerSize=2, Color="red");
            hold on;
            plot3([x1, x1 + scale*d12_x], [y1, y1 + scale*d12_y], [z1, z1 + scale*d12_z], Marker=">", LineStyle='-', MarkerSize=2, Color="red");
            hold on;
            plot3([x1, x1 + scale*d13_x], [y1, y1 + scale*d13_y], [z1, z1 + scale*d13_z], Marker=">", LineStyle='-', MarkerSize=2,  Color="red");
            hold on;
            % set axis lim. 
            grid on;
            xlim([x_lim1, x_lim2]); 
            ylim([y_lim1, y_lim2]); 
            zlim([z_lim1, z_lim2]);
            xlabel("x-axis"); 
            ylabel("y-axis"); 
            zlabel("z-axis")
            % title of plot. 
            title([title_str]);
            % hold on; 
        end

    end

end


