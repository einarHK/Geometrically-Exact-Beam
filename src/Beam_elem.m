
% single beam element class. 
classdef Beam_elem < handle
    properties
        x1 % initial coordinate - node 1. 
        x2 % initial coordinate - node 2. 
        L0 % initial beam length. 
        dof % degrees of freedom per node. 
        n_constraints % number of constraints per node. 
        d1_A % initial director 1 (shear, z) - node 1. 
        d2_A % initial director 2 (shear, y) - node 1. 
        d3_A % initial director 3 (axial) - node 1. 
        d1_B % initial director 1 (shear, z) - node 2. 
        d2_B % initial director 2 (shear, y) - node 2. 
        d3_B % initial director 3 (axial) - node 2. 
        lambda_At % lagrange multipliers - node 1. 
        lambda_Bt % lagrange multipliers - node 2. 
        x1_t % time dependent x1 pos - node 1. 
        x2_t % time dependent x2 pos - node 2. 
        v1_t % time dependent velocity vector - node 1. 
        v2_t % time dependent velocity vector - node 2. 
        d1_At % time dependent d1 director node 1. 
        d2_At % time dependent d2 director node 1. 
        d3_At % time dependent d3 director node 1
        d1_Bt % time dependent d1 director node 2. 
        d2_Bt % time dependent d2 director node 2. 
        d3_Bt % time dependent d3 director node 2. 
        w1_At % director velocity 1 node 1. 
        w2_At % director velocity 2 node 1. 
        w3_At % director velocity 3 node 1. 
        w1_Bt % director velocity 1 node 2. 
        w2_Bt % director velocity 2 node 2. 
        w3_Bt % director velocity 3 node 2. 
        fixed_dof % fixed degrees of freedom. 
        gamma_ref % reference axial strain vector. 
        omega_ref % reference curvature strain vector. 
        e0 % reference strain variable. 
        s0 % reference stress variable. 
    end

    methods
        % constructor.
        function obj = Beam_elem(x1, x2, L0, d3_A, d2_A, d1_A, d3_B, d2_B, d1_B, lambda_A, lambda_B, dof, n_constraints, fixed_dof, gamma_ref, omega_ref)
            obj.x1 = x1; 
            obj.x1_t = x1; 
            obj.x2 = x2; 
            obj.x2_t = x2;
            obj.L0 = L0;
            obj.d1_A = d1_A;
            obj.d2_A = d2_A;
            obj.d3_A = d3_A;
            obj.d1_At = d1_A; 
            obj.d2_At = d2_A;
            obj.d3_At = d3_A; 
            obj.d1_B = d1_B;
            obj.d2_B = d2_B; 
            obj.d3_B = d3_B;
            obj.d1_Bt = d1_B;
            obj.d2_Bt = d2_B;
            obj.d3_Bt = d3_B;
            obj.lambda_At = lambda_A;
            obj.lambda_Bt = lambda_B; 
            obj.dof = dof; 
            obj.n_constraints = n_constraints;
            obj.fixed_dof = fixed_dof;
            obj.gamma_ref = gamma_ref; 
            obj.omega_ref = omega_ref; 
            obj.v1_t = zeros(3,1); % assume 0 at start. 
            obj.v2_t = zeros(3,1); 
            obj.w1_At = zeros(3,1); 
            obj.w2_At = zeros(3,1); 
            obj.w3_At = zeros(3,1); 
            obj.w1_Bt = zeros(3,1); 
            obj.w2_Bt = zeros(3,1); 
            obj.w3_Bt = zeros(3,1); 
        end
        
        % returns Gauss quadrature points and weights. 
        function [weights, points] = gauss_quadrature(obj, n_points)
            if (n_points == 1)
                weights = [2]; 
                points = [0];
            elseif (n_points == 2)
                weights = [1, 1];
                points = [-1/sqrt(3), 1/sqrt(3)]; 
            else
                weights = [5/9, 8/9, 5/9]; 
                points = [-sqrt(0.6), 0, sqrt(0.6)];
            end
        end

        % returns X, Y, Z arrays containing the node position values.
       function [X, Y, Z] = compute_coordinates(obj)
            x_1 = obj.x1_t(1); % current x coordinate of node 1. 
            x_2 = obj.x2_t(1); % current x coordinate of node 2. 
            y1 = obj.x1_t(2); % current y coordinate of node 1. 
            y2 = obj.x2_t(2); % current y coordinate of node 2. 
            z1 = obj.x1_t(3); % current z coordinate of node 1. 
            z2 = obj.x2_t(3); % current z coordinate of node 2. 
            X = [x_1, x_2]; % x coordinates, node 1 and node 2. 
            Y = [y1, y2]; % y coordinates, node 1 and node 2. 
            Z = [z1, z2]; % z coordinates, node 1 and node 2. 
       end

       % computes the free degrees of freedom for the beam. 
       function free_dof = compute_free_dof(obj)
            free_dof = 1:2*obj.dof; 
            free_dof = setdiff(free_dof, obj.fixed_dof); % removes indices of fixed dof. 
       end

        % computes N1 = 0.5 * (1 - s) shape function. 
        function N1 = compute_N1(obj, s)
            N1 = 0.5 * (1 - s);
        end

        % computes N2 = 0.5 * (1 + s) shape function. 
        function N2 = compute_N2(obj, s)
            N2 = 0.5 * (1 + s);
        end
        
        % computes the derivative of the N1 shape function, w.r.t beam axis.
        function dN1 = compute_dN1(obj)
            dN1 = -(1/obj.L0);
        end

        % computes the derivative of the N2 shape function, w.r.t beam
        % axis.
        function dN2 = compute_dN2(obj)
            dN2 = (1/obj.L0); 
        end
        
        % computes the nodal value vector at first node. 
        function q1 = compute_q1(obj)
            q1 = [obj.x1_t;obj.d1_At;obj.d2_At;obj.d3_At];
        end
        
        % computes the nodal value vector at second node. 
        function q2 = compute_q2(obj)
            q2 = [obj.x2_t; obj.d1_Bt; obj.d2_Bt; obj.d3_Bt];
        end

        % computes q. 
        function q = compute_q(obj, s)
            q1 = obj.compute_q1();
            q2 = obj.compute_q2();
            n1 = obj.compute_N1(s);
            n2 = obj.compute_N2(s);
            I3 = eye(3);
            N1 = I3 * n1; 
            N2 = I3 * n2; 
            Z3 = zeros(3,3);
            Q = [N1, Z3, Z3, Z3, N2, Z3, Z3, Z3; Z3, N1, Z3, Z3, Z3, N2, Z3, Z3; Z3, Z3, N1, Z3, Z3, Z3, N2, Z3; Z3, Z3, Z3, N1, Z3, Z3, Z3, N2];
            q = Q * [q1;q2];            
        end

        % interpolation matrix. 
        % Assumes input vector of size 24x1. 
        function Q = compute_Q(obj, s)
            I3 = eye(3); 
            n1 = obj.compute_N1(s);
            n2 = obj.compute_N2(s);
            N1 = n1 * I3; 
            N2 = n2 * I3; 
            Z3 = zeros(3,3);
            Q1 = [N1, Z3, Z3, Z3, N2, Z3, Z3, Z3]; 
            Q2 = [Z3, N1, Z3, Z3, Z3, N2, Z3, Z3]; 
            Q3 = [Z3, Z3, N1, Z3, Z3, Z3, N2, Z3]; 
            Q4 = [Z3, Z3, Z3, N1, Z3, Z3, Z3, N2]; 
            Q = [Q1; Q2; Q3; Q4];
        end
        
        % differential interpolation operator. 
        % assumes input vector of size 6x1. 
        function dQ = compute_dQ(obj)
            dn1 = obj.compute_dN1();
            dn2 = obj.compute_dN2();
            I3 = eye(3); 
            Z3 = zeros(3,3);
            dN1 = dn1 * I3; 
            dN2 = dn2 * I3; 
            dQ1 = [dN1, Z3, Z3, Z3, dN2, Z3, Z3, Z3];
            dQ2 = [Z3, dN1, Z3, Z3, Z3, dN2, Z3, Z3]; 
            dQ3 = [Z3, Z3, dN1, Z3, Z3, Z3, dN2, Z3];
            dQ4 = [Z3, Z3, Z3, dN1, Z3, Z3, Z3, dN2];
            dQ = [dQ1; dQ2; dQ3; dQ4];
        end

        % computes gamma reference, when directors not aligned with
        % coordinate axes. 
        function gamma_ref = compute_gamma_ref(obj)
            dx = obj.x2 - obj.x1; 
            dx_size = sqrt(transpose(dx) * dx); 
            % beam axis unit vector. 
            gamma_ref = dx / dx_size; 
        end
        
        % computes the gamma vector. 
        % interpolation value s - in [-1,1]. 
        function gamma = compute_gamma(obj, s)
            n1 = obj.compute_N1(s); % n1 shape function. 
            n2 = obj.compute_N2(s); % n2 shape function. 
            dn1 = obj.compute_dN1(); % dn1/ds value. 
            dn2 = obj.compute_dN2(); % dn2/ds value. 
            dx = dn1 * obj.x1_t + dn2 * obj.x2_t; % the derivative of x w.r.t beam axis. 
            d1 = n1 * obj.d1_At + n2 * obj.d1_Bt; % interpolated d1 vector. 
            d2 = n1 * obj.d2_At + n2 * obj.d2_Bt; % interpolated d2 vector. 
            d3 = n1 * obj.d3_At + n2 * obj.d3_Bt; % interpolated d3 vector. 
            gamma1 = transpose(d1) * dx; % gamma 1 value.  
            gamma2 = transpose(d2) * dx; % gamma 2 value. 
            gamma3 = transpose(d3) * dx; % gamma 3 value. 
            gamma = [gamma1; gamma2; gamma3]; % gamma vector. 
        end
        
        % computes the omega vector. 
        % interpolation value s element of [-1,1].
        function omega = compute_omega(obj, s)
            n1 = obj.compute_N1(s); 
            n2 = obj.compute_N2(s);
            dn1 = obj.compute_dN1(); 
            dn2 = obj.compute_dN2();
            d1 = n1 * obj.d1_At + n2 * obj.d1_Bt; % interpolate the d_i vectors. 
            d2 = n1 * obj.d2_At + n2 * obj.d2_Bt;
            d3 = n1 * obj.d3_At + n2 * obj.d3_Bt;
            delta_d1 = dn1 * obj.d1_At + dn2 * obj.d1_Bt;
            delta_d2 = dn1 * obj.d2_At + dn2 * obj.d2_Bt;
            delta_d3 = dn1 * obj.d3_At + dn2 * obj.d3_Bt;
            omega1 = 0.5 * (transpose(d3)*delta_d2 - transpose(d2)*delta_d3);
            omega2 = 0.5 * (transpose(d1)*delta_d3 - transpose(d3)*delta_d1);
            omega3 = 0.5 * (transpose(d2)*delta_d1 - transpose(d1)*delta_d2);
            omega = [omega1; omega2; omega3]; 
        end
        
        % computes the strain vector e. 
        function e = compute_e(obj, s, gamma_ref, omega_ref)
            gamma = obj.compute_gamma(s); 
            omega = obj.compute_omega(s); 
            e = [gamma - gamma_ref; omega - omega_ref];
        end

        % computes the stresses for the given beam. 
        function s = compute_s(obj, s, C, gamma_ref, omega_ref)
            e = obj.compute_e(s, gamma_ref, omega_ref); 
            s = C * e; % linear constitutive material law. 
        end

        % computes B1 part of B matrix. 
        function B1 = compute_B1(obj, s)
            n1 = obj.compute_N1(s);
            n2 = obj.compute_N2(s);
            d1 = n1 * obj.d1_At + n2 * obj.d1_Bt;
            d2 = n1 * obj.d2_At + n2 * obj.d2_Bt;
            d3 = n1 * obj.d3_At + n2 * obj.d3_Bt;
            D = obj.compute_dQ(); % differential operator. 
            Z13 = zeros(1,3);
            B1 = 0.5 * [2*transpose(d1), Z13, Z13, Z13; 2*transpose(d2), Z13, Z13, Z13; 2*transpose(d3), Z13, Z13, Z13; Z13, Z13, transpose(d3), -transpose(d2); Z13, -transpose(d3), Z13, transpose(d1); Z13, transpose(d2), -transpose(d1), Z13]; 
            B1 = B1 * D; 
        end

        % computes B2 part of B matrix. 
        function B2 = compute_B2(obj, s)
            Z13 = zeros(1,3); 
            dn1 = obj.compute_dN1();
            dn2 = obj.compute_dN2();
            dx = dn1 * obj.x1_t + dn2 * obj.x2_t;
            B21 = [Z13, 2*transpose(dx), Z13, Z13]; 
            B22 = [Z13, Z13, 2*transpose(dx), Z13];
            B23 = [Z13, Z13, Z13, 2*transpose(dx)];
            delta_d1 = dn1 * obj.d1_At + dn2 * obj.d1_Bt;
            delta_d2 = dn1 * obj.d2_At + dn2 * obj.d2_Bt;
            delta_d3 = dn1 * obj.d3_At + dn2 * obj.d3_Bt; 
            B24 = [Z13, Z13, -transpose(delta_d3), transpose(delta_d2)]; 
            B25 = [Z13, transpose(delta_d3), Z13, -transpose(delta_d1)];
            B26 = [Z13, -transpose(delta_d2), transpose(delta_d1), Z13]; 
            B2 = 0.5 * [B21; B22; B23; B24; B25; B26]; 
            Q = obj.compute_Q(s); % calculate interpolation matrix. 
            B2 = B2 * Q; % multiply by interpolation matrix, such that takes input vector of 1x24. 
        end

        % computes the strain matrix B. 
        function B = compute_B(obj, s)
            B1 = obj.compute_B1(s); 
            B2 = obj.compute_B2(s); 
            B = B1 + B2; 
        end       

        % computes D_gamma differential operator. 
        % Input: 24x1 vector. 
        function D_gamma = compute_D_gamma(obj, s)
            dn1 = obj.compute_dN1();
            dn2 = obj.compute_dN2();
            I3 = eye(3); % 3x3 identity matrix.  
            dN1 = dn1 * I3; 
            dN2 = dn2 * I3; 
            N1 = obj.compute_N1(s); 
            N2 = obj.compute_N2(s);
            Z3 = zeros(3,3); 
            D11 = [dN1, Z3, Z3, Z3, dN2, Z3, Z3, Z3]; 
            D12 = [Z3, N1 * I3, Z3, Z3, Z3, N2 * I3, Z3, Z3]; 
            D13 = [Z3, Z3, N1 * I3, Z3, Z3, Z3, N2 * I3, Z3]; 
            D14 = [Z3, Z3, Z3, N1 * I3, Z3, Z3, Z3, N2 * I3]; 
            D_gamma = [D11; D12; D13; D14];
        end
        
        % computes U2_bar matrix. 
        function U2_bar = compute_U2_bar(obj, s_vector)
            s1 = s_vector(1);
            s2 = s_vector(2);
            s3 = s_vector(3); 
            I3 = eye(3); % 3x3 identity matrix. 
            Z3 = zeros(3,3); % 3x3 zero matrix. 
            U21 = [Z3, s1 * I3, s2 * I3, s3 * I3]; 
            U22 = [s1 * I3, Z3, Z3, Z3]; 
            U23 = [s2 * I3, Z3, Z3, Z3]; 
            U24 = [s3 * I3, Z3, Z3, Z3]; 
            U2_bar = [U21; U22; U23; U24]; 
        end

        % computes U2_hat. 
        function U2_hat = compute_U2_hat(obj, s_vector)
            s4 = s_vector(4); 
            s5 = s_vector(5); 
            s6 = s_vector(6); 
            Z3 = zeros(3,3); 
            I3 = eye(3); 
            U2_1 = [Z3, Z3, Z3, Z3, Z3, Z3, Z3, Z3]; 
            U2_2 = [Z3, Z3, Z3, Z3, Z3, -s6*I3, Z3, s5*I3]; 
            U2_3 = [Z3, Z3, Z3, s6*I3, Z3, Z3, -s5*I3, Z3];
            U2_4 = [Z3, Z3, s6*I3, Z3, Z3, Z3, Z3, -s4*I3];
            U2_5 = [Z3, Z3, Z3, Z3, Z3, Z3, Z3, Z3];
            U2_6 = [Z3, -s6*I3, Z3, Z3, Z3, Z3, s4*I3, Z3]; 
            U2_7 = [Z3, Z3, -s5*I3, Z3, Z3, s4*I3, Z3, Z3]; 
            U2_8 = [Z3, s5*I3, Z3, -s4*I3, Z3, Z3, Z3, Z3]; 
            U2_hat = 0.5 * [U2_1; U2_2; U2_3; U2_4; U2_5; U2_6;U2_7;U2_8];
        end

        function U2_hat_2 = compute_U2_hat_2(obj, s_vector)
            s4 = s_vector(4); 
            s5 = s_vector(5); 
            s6 = s_vector(6); 
            I3 = eye(3); 
            Z3 = zeros(3,3);
            U2_1 = [Z3, Z3, Z3, Z3, Z3, Z3, Z3, Z3]; 
            U2_2 = [Z3, Z3, Z3, Z3, Z3, Z3, -s6*I3, s5*I3]; 
            U2_3 = [Z3, Z3, Z3, Z3, Z3, s6*I3, Z3, -s4*I3];
            U2_4 = [Z3, Z3, Z3, Z3, Z3, -s5*I3, s4*I3, Z3]; 
            U2_5 = [Z3, Z3, Z3, Z3, Z3, Z3, Z3, Z3];
            U2_6 = [Z3, Z3, s6*I3, -s5*I3, Z3, Z3, Z3, Z3]; 
            U2_7 = [Z3, -s6*I3, Z3, s4*I3, Z3, Z3, Z3, Z3]; 
            U2_8 = [Z3, s5*I3, -s4*I3, Z3, Z3, Z3, Z3, Z3]; 
            U2_hat_2 = 0.5 * [U2_1; U2_2; U2_3; U2_4; U2_5; U2_6; U2_7; U2_8];
        end
        
        % computes D_omega. 
        function D_omega = compute_D_omega(obj, s)
            I12 = eye(12); 
            n1 = obj.compute_N1(s); 
            n2 = obj.compute_N2(s); 
            N1 = n1 * I12; 
            N2 = n2 * I12;
            dQ = obj.compute_dQ();
            D_omega = [N1, N2; dQ]; 
        end
        
        % computes U2(s) matrix - using gauss quadrature rule. 
        function U2 = compute_U2_integral(obj, n_gauss_points, C, gamma_ref, omega_ref) 
            [weights, points ] = obj.gauss_quadrature(n_gauss_points);
            U2 = zeros(2*obj.dof, 2*obj.dof);
            dx = obj.L0/2; 
            for i=1:n_gauss_points
                weight = weights(i); 
                point = points(i); 
                D_gamma = obj.compute_D_gamma(point);
                D_omega = obj.compute_D_omega(point); 
                s_vector = obj.compute_s(point, C, gamma_ref, omega_ref); 
                U2_bar = obj.compute_U2_bar(s_vector); 
                U2_hat = obj.compute_U2_hat(s_vector); 
                U2 = U2 + (transpose(D_gamma) * U2_bar * D_gamma + transpose(D_omega) * U2_hat * D_omega) * (dx * weight); 
            end
        end
        
        % derived U2 matrix. 
        function U2 = compute_U2_integral_2(obj, n_gauss_points, C, gamma_ref, omega_ref) 
            [weights, points ] = obj.gauss_quadrature(n_gauss_points);
            U2 = zeros(2*obj.dof, 2*obj.dof);
            dx = obj.L0/2; 
            for i=1:n_gauss_points
                weight = weights(i); 
                point = points(i); 
                D_gamma = obj.compute_D_gamma(point);
                D_omega = obj.compute_D_omega(point); 
                s_vector = obj.compute_s(point, C, gamma_ref, omega_ref); 
                U2_bar = obj.compute_U2_bar(s_vector); 
                U2_hat = obj.compute_U2_hat_2(s_vector); 
                U2 = (transpose(D_gamma) * U2_bar * D_gamma + transpose(D_omega) * U2_hat * D_omega) * (dx * weight); 
            end
        end
       
        
        % computes h - node 1. 
        function h1 = compute_h1(obj)
            d1 = obj.d1_At;
            d2 = obj.d2_At; 
            d3 = obj.d3_At; 
            h11 = transpose(d1) * d1 - 1; 
            h12 = transpose(d2) * d2 - 1; 
            h13 = transpose(d3) * d3 - 1; 
            h14 = 2 * transpose(d2) * d3;
            h15 = 2 * transpose(d1) * d3; 
            h16 = 2 * transpose(d1) * d2; 
            h1 = 0.5 * [h11; h12; h13; h14; h15; h16]; 
        end
        
        % computes h - node 2. 
        function h2 = compute_h2(obj)
            d1 = obj.d1_Bt; 
            d2 = obj.d2_Bt; 
            d3 = obj.d3_Bt; 
            h21 = transpose(d1) * d1 - 1; 
            h22 = transpose(d2) * d2 - 1; 
            h23 = transpose(d3) * d3 - 1; 
            h24 = 2 * transpose(d2) * d3; 
            h25 = 2 * transpose(d1) * d3; 
            h26 = 2 * transpose(d1) * d2; 
            h2 = 0.5 * [h21; h22; h23; h24; h25; h26];
        end
        
        % computes jacobian of the h vector - for node 1. 
        % computed at nodal level. 
        function H = compute_H_node1(obj)
            d1 = obj.d1_At; 
            d2 = obj.d2_At;
            d3 = obj.d3_At; 
            z13 = zeros(1,3); 
            H1 = [z13, transpose(d1), z13, z13];
            H2 = [z13, z13, transpose(d2), z13]; 
            H3 = [z13, z13, z13, transpose(d3)]; 
            H4 = [z13, z13, transpose(d3), transpose(d2)];
            H5 = [z13, transpose(d3), z13, transpose(d1)];
            H6 = [z13, transpose(d2), transpose(d1), z13];
            H = [H1; H2; H3; H4; H5; H6];
        end
        
        % computes the jacobian of the h vector - for node 2. 
        function H = compute_H_node2(obj)
            d1 = obj.d1_Bt;
            d2 = obj.d2_Bt; 
            d3 = obj.d3_Bt;
            z13 = zeros(1,3); 
            H1 = [z13, transpose(d1), z13, z13]; 
            H2 = [z13, z13, transpose(d2), z13];
            H3 = [z13, z13, z13, transpose(d3)];
            H4 = [z13, z13, transpose(d3), transpose(d2)];
            H5 = [z13, transpose(d3), z13, transpose(d1)];
            H6 = [z13, transpose(d2), transpose(d1), z13];
            H = [H1; H2; H3; H4; H5; H6]; 
        end

        % computes the total jacobian for the entire beam element. 
        function H = compute_H_tot(obj)
            H1 = obj.compute_H_node1(); % H matrix - node 1. 
            H2 = obj.compute_H_node2(); % H matrix - node 2. 
            [n_rows, n_cols] = size(H1); % [6, 12] 
            Z = zeros(n_rows, n_cols); % zero matrix 6x12. 
            H = [H1, Z; Z, H2]; % size 12 x 24
        end

        % computes the internal force vector. 
        function f_int = compute_f_int(obj, n_gauss_points, C, gamma_ref, omega_ref)
            [weights, points] = obj.gauss_quadrature(n_gauss_points);
            f_int = zeros(2*obj.dof, 1); % 24 x 1. 
            % f_int = zeros(2*obj.dof, 1); % initialize force vector - dimensions 12x1, 6x1 for each node. 
            dx = obj.L0/2; 
            for i=1:n_gauss_points
                weight = weights(i);
                point = points(i);
                % compute the stresses coming from this point. 
                s = obj.compute_s(point, C, gamma_ref, omega_ref); 
                % compute the strain matrix. 
                B = obj.compute_B(point); 
                % copute force contribution. 
                f_contrib = (transpose(B) * s) * (weight*dx); 
                % add to the total f_int. 
                f_int = f_int + f_contrib;
            end
        end
            
        % computes the force coming from the constraint at the discrete
        % nodes - node 1. 
        function f_H = compute_f_H1(obj)
            H = obj.compute_H_node1(); % compute the jacobian of h at node 1.
            f_H = transpose(H) * obj.lambda_At; % compute the force. 
        end

        % computes the force coming from the constraint at the discrete
        % nodes - node 2.
        function f_H = compute_f_H2(obj)
            H = obj.compute_H_node2(); % compute the jacobian of h at node 1.
            f_H = transpose(H) * obj.lambda_Bt; % compute the force. 
        end 

        % computes the V matrix - for node 1. 
        function V1 = compute_V1(obj)
            lambda = obj.lambda_At; % lambda vector at node 1. 
            v1 = lambda(1); 
            v2 = lambda(2); 
            v3 = lambda(3); 
            v4 = lambda(4); 
            v5 = lambda(5); 
            v6 = lambda(6);
            Z3 = zeros(3,3); 
            I3 = eye(3);
            V11 = [Z3, Z3, Z3, Z3];
            V12 = [Z3, v1 * I3, v6 * I3, v5 * I3]; 
            V13 = [Z3, v6 * I3, v2 * I3, v4 * I3]; 
            V14 = [Z3, v5 * I3, v4 * I3, v3 * I3];
            V1 = [V11; V12; V13; V14]; 
        end

        
        % computes the (dH^T/dq * v) matrix for node 2. 
        function V2 = compute_V2(obj)
            Z3 = zeros(3,3); 
            lambda = obj.lambda_Bt; 
            v1 = lambda(1);
            v2 = lambda(2); 
            v3 = lambda(3); 
            v4 = lambda(4);
            v5 = lambda(5); 
            v6 = lambda(6); 
            I3 = eye(3);
            V21 = [Z3, Z3, Z3, Z3];
            V22 = [Z3, v1 * I3, v6 * I3, v5 * I3]; 
            V23 = [Z3, v6 * I3, v2 * I3, v4 * I3]; 
            V24 = [Z3, v5 * I3, v4 * I3, v3 * I3]; 
            V2 = [V21; V22; V23; V24]; 
        end

        % computes the total dH^T/dq matrix for the beam element. 
        function V = compute_V_tot(obj)
            V1 = obj.compute_V1(); % node 1. 
            V2 = obj.compute_V2(); % node 2. 
            Z12 = zeros(obj.dof, obj.dof); % zero matrix of size 12x12. 
            V = [V1, Z12; Z12, V2]; 
        end

        % interpolation matrix - lambda nodal variables. 
        function Q_lambda = compute_Q_lambda(obj, s)
            n1 = obj.compute_N1(s); 
            n2 = obj.compute_N2(s); 
            I6 = eye(6); 
            Q_lambda = [n1 * I6, n2 * I6];
        end
        
        % computes material tanget stiffness matrix. 
        function Kt_m = compute_Kt_m(obj, n_gauss_points, C)
            [weights, points] = obj.gauss_quadrature(n_gauss_points);
            Kt_m = zeros(2*obj.dof, 2*obj.dof); 
            dx = (obj.L0)/2; % differential of the integral domain. 
            for i=1:n_gauss_points
                weight = weights(i); % scalar weight. 
                point = points(i); % integration point. 
                B = obj.compute_B(point); % computes the strain matrix. 
                Kt_m = Kt_m + (transpose(B)*C*B) * (weight * dx); % add to the Kt_m matrix. 
            end
        end

        % computes the Kt stiffness tangent matrix, related to internal
        % force term. 
        function Kt = compute_Kt(obj, n_gauss_points, C, gamma_ref, omega_ref)
            % compute material tangent stiffness matrix. 
            Kt_m = obj.compute_Kt_m(n_gauss_points, C); 
            % compute geometric tangent stiffness matrix. 
            Kt_g = obj.compute_U2_integral_2(n_gauss_points, C, gamma_ref, omega_ref); % derived U2 matrix. 
            % Kt_g = obj.compute_U2_integral(n_gauss_points, C, gamma_ref, omega_ref);
            % add the matrices. 
            Kt = Kt_g + Kt_m;
        end

        % computes the S matrix. 
        % dimensions - 36 x 36. 
        function S = compute_S_mat(obj, n_gauss_points, C, gamma_ref, omega_ref)
            Kt = obj.compute_Kt(n_gauss_points, C, gamma_ref, omega_ref);
            dHv_dq = obj.compute_V_tot();  
            S11 = Kt + dHv_dq; % size 24x24
            S21 = obj.compute_H_tot(); % size 24x12
            S12 = transpose(S21); % size 12x24
            [rows_S12, cols_S12] = size(S12);
            [rows_S21, cols_S21] = size(S21); 
            S22 = zeros(rows_S21, cols_S12); % size 24x12.
            S = [S11, S12; S21, S22];
        end
        
        
        % assumes e and s have dimensions 6 for each element - not as nodal
        % point values. 
        function S = compute_full_S_mat(obj, s_point, n_gauss_points, C, gamma_ref, omega_ref)
            % r = B(q)^T * s + H(q)^T * X - f_ext.
            % derivative of r w.r.t q. 
            dr_dq = obj.compute_U2_integral_2(n_gauss_points, C, gamma_ref, omega_ref) + obj.compute_V_tot(); 
            % derivative of r w.r.t e.
            dr_de = zeros(2 * obj.dof, 6); % assuming e has dim 6. 
            % derivative of r w.r.t s. 
            dx = obj.L0 / 2; 
            [weights, points] = obj.gauss_quadrature(n_gauss_points);
            dr_ds = zeros(24,6); 
            for i=1:n_gauss_points
                point = points(i); 
                weight = weights(i); 
                B = obj.compute_B(point); 
                dr_ds = dr_ds + transpose(B) * (dx * weight);
            end
            % derivative of r w.r.t lambda. 
            H = obj.compute_H_tot();
            dr_dLambda = transpose(H);
            S1 = [dr_dq, dr_de, dr_ds, dr_dLambda]; 
            % derivative of f(e,e(q)) = e - e(q) w.r.t e. 
            df_de = eye(6); 
            % derivative of f w.r.t q. 
            df_dq = -obj.compute_B(s_point); 
            % derivative of f w.r.t s. 
            df_ds = zeros(6,6); 
            % derivative of f w.r.t lambda. 
            df_dLambda = zeros(6,2 * obj.n_constraints); 
            % S2 row. 
            S2 = [df_dq, df_de, df_ds, df_dLambda]; 
            % derivative of g(e,s) w.r.t q. 
            dg_dq = zeros(6, 2 * obj.dof);
            % derivative of g(e,s) w.r.t e. 
            dg_de = -C; 
            % derivative of g(e,s) w.r.t s. 
            dg_ds = eye(6); 
            % derivative of g(e,s) w.r.t lambda. 
            dg_dLambda = zeros(6, 2 * obj.n_constraints); 
            % S3 row. 
            S3 = [dg_dq, dg_de, dg_ds, dg_dLambda]; 
            % derivative of h w.r.t q. 
            dh_dq = obj.compute_H_tot(); 
            % derivative of h w.r.t e. 
            dh_de = zeros(2 * obj.n_constraints, 6); 
            % derivative of h w.r.t s. 
            dh_ds = zeros(2 * obj.n_constraints, 6); 
            % derivative of h w.r.t lambda. 
            dh_dLambda = zeros(2 * obj.n_constraints, 2*obj.n_constraints);
            % S4 row. 
            S4 = [dh_dq, dh_de, dh_ds, dh_dLambda]; 
            % S matrix. 
            S = [S1; S2; S3; S4]; 
        end

        % computes the residual vector for the forces. 
        % must be 24x1 vector. 
        function residual = compute_residual(obj, n_gauss_points, C, gamma_ref, omega_ref, f_ext)
            % compute internal force term. 
            f_int = obj.compute_f_int(n_gauss_points, C, gamma_ref, omega_ref); 
            % compute constraint force term - at each node. 
            f_H1 = obj.compute_f_H1(); 
            f_H2 = obj.compute_f_H2(); 
            f_H = [f_H1; f_H2]; 
            % compute the residual. 
            residual = f_int + f_H - f_ext;
        end

        % computes the h(q) value for the constraint vector for given nodal
        % values. dimensions 12x1, since each node has 6 constraints. 
        function h_q = compute_h_q(obj)
            h1 = obj.compute_h1(); % restriction - node 1. 
            h2 = obj.compute_h2(); % restriction - node 2.  
            h_q = [h1; h2]; 
        end
       
        % computes the right hand side vector for the linearized system of
        % equations. dimensions are (24 + 12) x 1 = 36x1. 
        function rhs = compute_rhs(obj, n_gauss_points, C, gamma_ref, omega_ref, f_ext)
            % residual due to force terms. 
            residual = obj.compute_residual(n_gauss_points, C, gamma_ref, omega_ref, f_ext); 
            % error due to constraint terms. 
            h_q = obj.compute_h_q(); % error coming from constraint equation. 
            % right hand side vector. 
            rhs = [-residual; -h_q];
        end

        % updates the parameters for the nodes. 
        function update_params(obj, dx1, dx2, delta_d1_A, delta_d2_A, delta_d3_A, delta_d1_B, delta_d2_B, delta_d3_B)
            obj.x1_t = obj.x1_t + dx1; 
            obj.x2_t = obj.x2_t + dx2; 
            obj.d1_At = obj.d1_At + delta_d1_A;
            obj.d2_At = obj.d2_At + delta_d2_A;
            obj.d3_At = obj.d3_At + delta_d3_A;
            obj.d1_Bt = obj.d1_Bt + delta_d1_B;
            obj.d2_Bt = obj.d2_Bt + delta_d2_B;
            obj.d3_Bt = obj.d3_Bt + delta_d3_B;
        end

        % updates lambda vector. 
        function update_lambda(obj, delta_lambda_A, delta_lambda_B)
            obj.lambda_At = obj.lambda_At + delta_lambda_A;
            obj.lambda_Bt = obj.lambda_Bt + delta_lambda_B;
        end
        
         % updates the values of the velocity parameters. 
        function update_velocity_params(obj, v1, v2, w1_At, w2_At, w3_At, w1_Bt, w2_Bt, w3_Bt)
            obj.v1_t = v1; % velocity at node 1. 
            obj.v2_t = v2; % velocity at node 2. 
            obj.w1_At = w1_At; % rotation velocity node 1, director 1. 
            obj.w2_At = w2_At; % rotation velocity node 1, director 2. 
            obj.w3_At = w3_At; % rotation velocity node 1, director 3. 
            obj.w1_Bt = w1_Bt; % rotation velocity node 2, director 1. 
            obj.w2_Bt = w2_Bt; % rotation velocity node 2, director 2. 
            obj.w3_Bt = w3_Bt; % rotation velocity node 2, director 3. 
        end

        % updates the configuration of the beam. 
        function update_config(obj, delta_u)
            dx1 = delta_u(1:3); 
            delta_d1_A = delta_u(4:6); 
            delta_d2_A = delta_u(7:9);
            delta_d3_A = delta_u(10:12); 
            dx2 = delta_u(13:15); 
            delta_d1_B = delta_u(16:18); 
            delta_d2_B = delta_u(19:21); 
            delta_d3_B = delta_u(22:24); 
            delta_lambda_A = delta_u(25:30); 
            delta_lambda_B = delta_u(31:36); 
            obj.update_params(dx1, dx2, delta_d1_A, delta_d2_A, delta_d3_A, delta_d1_B, delta_d2_B, delta_d3_B);
            obj.update_lambda(delta_lambda_A, delta_lambda_B);
        end

        % displays current node positions. 
        function display_node_pos(obj)
            x1 = obj.x1_t(1); 
            y1 = obj.x1_t(2); 
            z1 = obj.x1_t(3); 
            x2 = obj.x2_t(1); 
            y2 = obj.x2_t(2); 
            z2 = obj.x2_t(3); 
            output1 = "Node 1 coords: (" + num2str(x1) + "m " + num2str(y1) + "m " + num2str(z1) + "m)";
            output2 = "Node 2 coords: (" + num2str(x2) + "m " + num2str(y2) + "m " + num2str(z2) + "m)";
            disp([output1]);
            disp([output2]);
        end
        
        % displays node directors. 
        function display_node_directors(obj)
            d1_A = obj.d1_At;
            d2_A = obj.d2_At;
            d3_A = obj.d3_At;
            d1_B = obj.d1_Bt; 
            d2_B = obj.d2_Bt; 
            d3_B = obj.d3_Bt;
            output1 = "Node 1 directors";
            output2 = "d1 = [" + num2str(d1_A(1)) + ", " + num2str(d1_A(2)) + ", " + num2str(d1_A(3)) + "]";
            output3 = "d2 = [" + num2str(d2_A(1)) + ", " + num2str(d2_A(2)) + ", " + num2str(d2_A(3)) + "]";
            output4 = "d3 = [" + num2str(d3_A(1)) + ", " + num2str(d3_A(2)) + ", " + num2str(d3_A(3)) + "]";
            output5 = "Node 2 directors:";
            output6 = "d1 = [" + num2str(d1_B(1)) + ", " + num2str(d1_B(2)) + ", " + num2str(d1_B(3)) + "]"; 
            output7 = "d2 = [" + num2str(d2_B(1)) + ", " + num2str(d2_B(2)) + ", " + num2str(d2_B(3)) + "]"; 
            output8 = "d3 = [" + num2str(d3_B(1)) + ", " + num2str(d3_B(2)) + ", " + num2str(d3_B(3)) + "]"; 
            disp(output1);
            disp(output2);
            disp(output3); 
            disp(output4); 
            disp(output5); 
            disp(output6);
            disp(output7);
            disp(output8);
        end

        % plots the current beam configuration. 
        function show_config(obj, x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scaling, title_str, scale_factor)
            [X, Y, Z] = obj.compute_coordinates(); % returns node 1 and node 2 X,Y,Z vals. 
            X = round(X, 5); 
            Y = round(Y, 5); 
            Z = round(Z, 5);
            % plot beam axis. 
            plot3(X, Y, Z, LineStyle='-', Color="red", Marker="o", MarkerFaceColor="blue", MarkerSize=2, LineWidth=1);
            hold on; 
            % plot directors. 
            x1 = obj.x1_t(1); 
            y1 = obj.x1_t(2); 
            z1 = obj.x1_t(3); 
            x2 = obj.x2_t(1); 
            y2 = obj.x2_t(2); 
            z2 = obj.x2_t(3); 
            % directors at first node. 
            d11_x = obj.d1_At(1);
            d11_y = obj.d1_At(2); 
            d11_z = obj.d1_At(3); 
            d12_x = obj.d2_At(1); 
            d12_y = obj.d2_At(2); 
            d12_z = obj.d2_At(3);
            d13_x = obj.d3_At(1); 
            d13_y = obj.d3_At(2); 
            d13_z = obj.d3_At(3);
            % directors at second node. 
            d21_x = obj.d1_Bt(1);
            d21_y = obj.d1_Bt(2); 
            d21_z = obj.d1_Bt(3); 
            d22_x = obj.d2_Bt(1); 
            d22_y = obj.d2_Bt(2); 
            d22_z = obj.d2_Bt(3);
            d23_x = obj.d3_Bt(1); 
            d23_y = obj.d3_Bt(2); 
            d23_z = obj.d3_Bt(3);
            if (scaling ~= 0) 
                dx = x2 - x1; 
                scale = sqrt(dx' * dx); 
            else
                scale = 1; 
            end
            scale = scale_factor * scale; 
            % plot directors at each node. 
            plot3([x1, x1 + scale * d11_x], [y1, y1 + scale *d11_y], [z1, z1 + scale*d11_z], Marker=">", LineStyle='-', Color="blue", MarkerSize=2, LineWidth=1);
            hold on;
            plot3([x1, x1 + scale*d12_x], [y1, y1 + scale*d12_y], [z1, z1 + scale*d12_z], Marker=">", LineStyle='-', Color="blue", MarkerSize=2, LineWidth=1);
            hold on;
            plot3([x1, x1 + scale*d13_x], [y1, y1 + scale*d13_y], [z1, z1 + scale*d13_z], Marker=">", LineStyle='-', Color="blue", MarkerSize=2, LineWidth=1);
            grid on;
            plot3([x2, x2 + scale*d21_x], [y2, y2 + scale*d21_y], [z2, z2 + scale*d21_z], Marker=">", LineStyle='-', Color="blue", MarkerSize=2, LineWidth=1);
            hold on;
            plot3([x2, x2 + scale*d12_x], [y2, y2 + scale*d22_y], [z2, z2 + scale*d22_z], Marker=">", LineStyle='-', Color="blue", MarkerSize=2, LineWidth=1);
            hold on;
            plot3([x2, x2 + scale*d23_x], [y2, y2 + scale*d23_y], [z2, z2 + scale*d23_z], Marker=">", LineStyle='-', Color="blue", MarkerSize=2, LineWidth=1);
            hold on;
            grid on;
            xlabel("x-axis"); 
            ylabel("y-axis");
            zlabel("z-axis");
            % set axis lim. 
            xlim([x_lim1, x_lim2]); 
            ylim([y_lim1, y_lim2]); 
            zlim([z_lim1, z_lim2]);
            % title of plot. 
            title([title_str]);
        end

    end
end
            



