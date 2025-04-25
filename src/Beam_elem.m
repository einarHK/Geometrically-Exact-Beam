
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
        mu % the lagrange multipliers for the (e - e(q)) equation. dim(muh) = 6
        chi % the lagrange multipliers for the balance equation. dim(chi) = dim(q) = 12. 
        gamma % the lagrange multipliers for the restriction eq. h(q) = 0. dim(gamma) = 6.
        xi % the lagrange multipliers for the constitutive eq. g = 0. dim(xi) = 6. 
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
            % init the Lagrange multipliers. 
            obj.mu = zeros(6,1); % for the compatibility eq. e - e(q) = 0. 
            obj.chi = zeros(24,1); % for the balance eq. 
            obj.gamma = zeros(12,1); % for the eq. h(q) = 0. 
            obj.xi = zeros(6,1); % for the constitutive eq. g=0. 
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
        
        % computes the derivative of B1(q)*a w.r.t q. 
        % should be 6x24, a is 12x1, already interpolated.  
        function dB1_dq = compute_dB1_dq(obj, s, a)
            Z3 = zeros(1,3); 
            Q = obj.compute_Q(s); % interpolation matrix. 
            a = Q * a; % interpolate a. 
            a0 = transpose(a(1:3)); 
            a1 = transpose(a(4:6)); 
            a2 = transpose(a(7:9)); 
            a3 = transpose(a(10:12)); 
            dB1 = [Z3, a0, Z3, Z3]; 
            dB2 = [Z3, Z3, a0, Z3]; 
            dB3 = [Z3, Z3, Z3, a0]; 
            dB4 = [Z3, Z3, -0.5 * a3, 0.5 * a2]; 
            dB5 = [Z3, 0.5 * a3, Z3, -0.5 * a1]; 
            dB6 = [Z3, -0.5 * a2, 0.5 * a1, Z3]; 
            dQ = obj.compute_dQ(); 
            dB1_dq = [dB1; dB2; dB3; dB4; dB5; dB6] * dQ; 
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

        % computes the derivative of B2(q) * a w.r.t q. 
        % a is 24x1 vector, s is the interpolation value. 
        function dB2_dq = compute_dB2_dq(obj, s, a)
            Z3 = zeros(1,3); 
            dQ = obj.compute_dQ(); 
            % interpolation matrix for the input vector. 
            Q = obj.compute_Q(s); 
            da_dq = dQ * a; % calculate the derivative of a. 
            da0_dq = transpose(da_dq(1:3)); 
            da1_dq = transpose(da_dq(4:6)); 
            da2_dq = transpose(da_dq(7:9)); 
            da3_dq = transpose(da_dq(10:12)); 
            dB1 = [Z3, da0_dq, Z3, Z3]; 
            dB2 = [Z3, Z3, da0_dq, Z3]; 
            dB3 = [Z3, Z3, Z3, da0_dq]; 
            dB4 = [Z3, Z3, -0.5 * da3_dq, 0.5 * da2_dq]; 
            dB5 = [Z3, da3_dq, Z3, -da1_dq]; 
            dB6 = [Z3, -da2_dq, da1_dq, Z3]; 
            dB2_dq = [dB1; dB2; dB3; dB4; dB5; dB6] * Q; 
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

        % computes U1_bar matrix - component for computing the d/dq(B(q)a),
        % where a is a 12x1 vector - 24 x1 when considering linear
        % interpolation - here s val is the interpolation value. 
        % a is a 24x1 vector. 
        function U1_bar = compute_U1_bar(obj, a)
            Z13 = zeros(1, 3); 
            dQ = obj.compute_dQ(); 
            % compute the derivative of a. 
            delta_a = dQ * a; 
            delta_a0 = delta_a(1:3); 
            delta_a1 = delta_a(4:6); 
            delta_a2 = delta_a(7:9);
            delta_a3 = delta_a(10:12); 
            U11 = [Z13, transpose(delta_a0), Z13, Z13]; 
            U12 = [Z13, Z13, transpose(delta_a0), Z13]; 
            U13 = [Z13, Z13, Z13, transpose(delta_a0)]; 
            U14 = [Z13, Z13, -0.5 * transpose(delta_a3), transpose(delta_a2)]; 
            U15 = [Z13, 0.5 * transpose(delta_a3), Z13, -0.5 * transpose(delta_a1)]; 
            U16 = [Z13, -0.5 * transpose(delta_a2), 0.5 * transpose(delta_a1), Z13]; 
            U1_bar = [U11; U12; U13; U14; U15; U16]; 
        end
        
        % assume a is 24x1 vector. 
        function U1_hat = compute_U1_hat(obj, s_val, a)
            n1 = obj.compute_N1(s_val); 
            n2 = obj.compute_N2(s_val);
            Z13 = zeros(1,3); 
            N1 = eye(12) * n1; 
            N2 = eye(12) * n2; 
            N_mat = [N1, N2]; 
            a_interp = N_mat * a; % interpolated a vector. 
            a0 = a_interp(1:3); 
            a1 = a_interp(4:6); 
            a2 = a_interp(7:9); 
            a3 = a_interp(10:12); 
            U11 = [transpose(a1), Z13, Z13, Z13]; 
            U12 = [transpose(a2), Z13, Z13, Z13]; 
            U13 = [transpose(a3), Z13, Z13, Z13]; 
            U14 = [Z13, Z13, 0.5 * transpose(a3), -0.5 * transpose(a2)]; 
            U15 = [Z13, -0.5 * transpose(a3), Z13, 0.5 * transpose(a1)]; 
            U16 = [Z13, 0.5 * transpose(a2), -0.5 * transpose(a1), Z13]; 
            U1_hat = [U11; U12; U13; U14; U15; U16]; 
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
            U2_bar = [U21; U22; U23; U24]; % 12x12 matrix. 
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
        
        % derived U2 matrix - computes the geometric tangent stiffness matrix integral. 
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

        % returns the cross product matrix for a vector d.
        % assuming the order of the cross product with d is x d - left
        % sided. 
        function d_mat = compute_cross_prod_mat(obj, d)
            d1 = d(1); 
            d2 = d(2); 
            d3 = d(3); 
            d_mat = [0, d3, -d2; -d3, 0, d1; d2, -d1, 0]; 
        end

        % Null space basis - node 1. 
        function N1 = compute_N1_basis(obj)
            I = eye(3); 
            Z3 = zeros(3,3); 
            d1_mat = obj.compute_cross_prod_mat(obj.d1_At); 
            d2_mat = obj.compute_cross_prod_mat(obj.d2_At); 
            d3_mat = obj.compute_cross_prod_mat(obj.d3_At); 
            N1 = [I, Z3, Z3, Z3; Z3, d1_mat, d2_mat, d3_mat]; 
            N1 = transpose(N1);
        end

        % Null space basis - node 2. 
        function N2 = compute_N2_basis(obj)
            I = eye(3); 
            Z3 = zeros(3,3); 
            d1_mat = obj.compute_cross_prod_mat(obj.d1_Bt); 
            d2_mat = obj.compute_cross_prod_mat(obj.d2_Bt); 
            d3_mat = obj.compute_cross_prod_mat(obj.d3_Bt); 
            N2 = [I, Z3, Z3, Z3; Z3, d1_mat, d2_mat, d3_mat]; 
            N2 = transpose(N2); 
        end

        % Compues the full null basis matrix. 
        function N = compute_N_basis(obj)
            N1 = obj.compute_N1_basis(); 
            N2 = obj.compute_N2_basis(); 
            N = [N1; N2]; 
        end

        % Computes the W1(q,a) = d(N(q)^T * a)/dq 
        % assume a is 1x24
        function W1 = compute_W1(obj, a)
            a1 = a(1:3); 
            a2 = a(4:6); 
            a3 = a(7:9); 
            a4 = a(10:12); 
            a5 = a(13:15); 
            a6 = a(16:18); 
            a7 = a(19:21); 
            a8 = a(22:24); 
            Z3 = zeros(3,3); 
            a1_mat = obj.compute_cross_prod_mat(a1); 
            a2_mat = obj.compute_cross_prod_mat(a2); 
            a3_mat = obj.compute_cross_prod_mat(a3); 
            a4_mat = obj.compute_cross_prod_mat(a4); 
            a5_mat = obj.compute_cross_prod_mat(a5); 
            a6_mat = obj.compute_cross_prod_mat(a6); 
            a7_mat = obj.compute_cross_prod_mat(a7); 
            a8_mat = obj.compute_cross_prod_mat(a8);
            W1 = -[Z3, Z3, Z3, Z3, Z3, Z3, Z3, Z3; Z3, a2_mat, a3_mat, a4_mat, Z3, a6_mat, a7_mat, a8_mat]; 
        end
        
        % computes W2 = d(N(q)^T * b)/dq, b - 1x6 vector. 
        function W2 = compute_W2(obj, b)
            b1 = b(1:3); 
            b2 = b(4:6); 
            b2_mat = obj.compute_cross_prod_mat(b2); 
            Z3 = zeros(3,3); 
            W2 = [Z3, Z3, Z3, Z3; Z3, b2_mat, Z3, Z3; Z3, Z3, b_mat, Z3; Z3, Z3, Z3, b2_mat]; 
        end

        % computes S_fix - using null basis. 
        function S = compute_S_fix_nullbasis(obj, n_gauss_points, C)
            % compute S_qq. 
            d_h_gamma_dq = obj.compute_d_H_gamma_dq();
            % U2(q,s) = d(B(q)^T * s)/dq
            U2 = obj.compute_U2_integral_2(n_gauss_points, C, obj.gamma_ref, obj.omega_ref); 
            % - d(B(q)^T * mu). 
            dB_mu_dq = -1 * obj.compute_dB_mu_dq(0); 
            % W2 = d(N(q) * chi)/dq
            
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
        
        function d_H_gammma_dq = compute_d_H_gamma_dq(obj)
            Z12 = zeros(12, 12); % 12 by 12 zero matrix. 
            Z3 = zeros(3,3); 
            I3 = eye(3); 
            g1 = obj.gamma(1); 
            g2 = obj.gamma(2); 
            g3 = obj.gamma(3); 
            g4 = obj.gamma(4); 
            g5 = obj.gamma(5); 
            g6 = obj.gamma(6); 
            g11 = obj.gamma(7); 
            g22 = obj.gamma(8); 
            g33 = obj.gamma(9); 
            g44 = obj.gamma(10); 
            g55 = obj.gamma(11); 
            g66 = obj.gamma(12);
            dH_gamma_dq_1 = [Z3, Z3, Z3, Z3; Z3, g1*I3, g6*I3, g5*I3; Z3, g6*I3, g2*I3, g4*I3; Z3, g5*I3, g4*I3, g3*I3]; 
            dH_gamma_dq_2 = [Z3, Z3, Z3, Z3; Z3, g11*I3, g66*I3, g55*I3; Z3, g66*I3, g22*I3, g44*I3; Z3, g55*I3, g44*I3, g33*I3];
            d_H_gammma_dq = [dH_gamma_dq_1, Z12; Z12, dH_gamma_dq_2]; 
        end

        % S components for the KKT matrix.
        % S_e_tilde_e_tilde matrix. 
        function S_e_tilde_e_tilde = compute_S_e_tilde_e_tilde(obj, C)
            S_e_tilde_e_tilde = transpose(C); 
        end

        function S_e_tilde_s_tilde = compute_S_e_tilde_s_tilde(obj, C)
            S_e_tilde_s_tilde = zeros(6,6); 
        end
        
        function S_e_tilde_x = compute_S_e_tilde_x(obj, C)
            S_e_tilde_q = zeros(6, 24); 
            S_e_tilde_e = -transpose(C); 
            S_e_tilde_s = zeros(6,6); 
            S_e_tilde_lambda = zeros(6, 12); 
            S_e_tilde_chi = zeros(6, 24); 
            S_e_tilde_mu = zeros(6,6); 
            S_e_tilde_gamma = zeros(6, 12); 
            S_e_tilde_x = [S_e_tilde_q, S_e_tilde_e, S_e_tilde_s, S_e_tilde_lambda, S_e_tilde_chi, S_e_tilde_mu, S_e_tilde_gamma];
        end

        function S_e_tilde_xi = compute_S_e_tilde_xi(obj, C, e_tilde, s_tilde)
            dg_de_tilde = obj.compute_dg_de_tilde(C, e_tilde, s_tilde); 
            S_e_tilde_xi = transpose(dg_de_tilde); 
        end

        function S_s_tilde_e_tilde = compute_S_s_tilde_e_tilde(obj)
            S_s_tilde_e_tilde = zeros(6,6); 
        end
        
        function S_s_tilde_s_tilde = compute_S_s_tilde_s_tilde(obj, C)
            S_s_tilde_s_tilde = transpose(C); 
        end

        function S_s_tilde_x = compute_S_s_tilde_x(obj, C)
            S_s_tilde_q = zeros(6, 24); 
            S_s_tilde_e = zeros(6, 6); 
            S_s_tilde_s = -transpose(C); 
            S_s_tilde_lambda = zeros(6, 12); 
            S_s_tilde_chi = zeros(6,24); 
            S_s_tilde_mu = zeros(6,6);
            S_s_tilde_gamma = zeros(6,12); 
            S_s_tilde_x = [S_s_tilde_q, S_s_tilde_e, S_s_tilde_s, S_s_tilde_lambda, S_s_tilde_chi, S_s_tilde_mu, S_s_tilde_gamma]; 
        end

        function S_s_tilde_xi = compute_S_s_tilde_xi(obj, C, e_tilde, s_tilde)
            dg_ds_tilde = obj.compute_dg_de_tilde(C, e_tilde, s_tilde); 
            S_s_tilde_xi = transpose(dg_ds_tilde); 
        end
        
        % calculates the derivative of B^T * mu, w.r.t q. 
        function dB_mu_dq = compute_dB_mu_dq(obj, s)
            U2_bar = obj.compute_U2_bar(obj.mu); 
            U2_hat = obj.compute_U2_hat_2(obj.mu); 
            D_gamma = obj.compute_D_gamma(s); 
            D_omega = obj.compute_D_omega(s); 
            dB_mu_dq = transpose(D_gamma) * U2_bar * D_gamma + transpose(D_omega) * U2_hat * D_omega;
        end

        function S_qq = compute_S_qq(obj)
            % compute d(H(q)^T * gamma)/dq
            dH_gamma_dq = obj.compute_d_H_gamma_dq();
            % compute -d(B(q)^T * mu)/dq - s = 0 interpolation value. 
            dB_mu_dq = obj.compute_dB_mu_dq(0); 
            % add the two matrices. 
            S_qq = dH_gamma_dq - dB_mu_dq; %  
        end
        
        function S_qe = compute_S_qe(obj)
            S_qe = zeros(24, 6); % 24x6 zero matrix. 
        end
        
        function dU_bar_dv = compute_dU_bar_dv(obj, v)
            v1 = v(1); 
            v2 = v(2); 
            v3 = v(3); 
            v4 = v(4); 
            v5 = v(5); 
            v6 = v(6); 
            v7 = v(7); 
            v8 = v(8); 
            v9 = v(9); 
            v10 = v(10);
            v11 = v(11); 
            v12 = v(12); 
            U1 = [v4, v7, v10, 0, 0, 0]; 
            U2 = [v5, v8, v11, 0, 0, 0]; 
            U3 = [v6, v9, v12, 0, 0, 0]; 
            U4 = [v1, 0, 0, 0, 0, 0]; 
            U5 = [v2, 0, 0, 0, 0, 0]; 
            U6 = [v3, 0, 0, 0, 0, 0]; 
            U7 = [0, v1, 0, 0, 0, 0]; 
            U8 = [0, v2, 0, 0, 0, 0]; 
            U9 = [0, v3, 0, 0, 0, 0]; 
            U10 = [0, 0, v1, 0, 0, 0]; 
            U11 = [0, 0, v2, 0, 0, 0]; 
            U12 = [0, 0, v3, 0, 0, 0]; 
            dU_bar_dv = [U1; U2; U3; U4; U5; U6; U7; U8; U9; U10; U11; U12];
        end
    
        function dU_hat_dv = compute_dU_hat_dv(obj, v) 
            U1_3 = zeros(3, 6);  
            U4 = [0, 0, 0, 0, v(22), -v(19)]; 
            U5 = [0, 0, 0, 0, v(23), -v(20)]; 
            U6 = [0, 0, 0, 0, v(24), -v(21)]; 
            U7 = [0, 0, 0, -v(22), 0, v(16)]; 
            U8 = [0, 0, 0, -v(23), 0, v(17)]; 
            U9 = [0, 0, 0, -v(24), 0, v(18)]; 
            U10 = [0, 0, 0, v(19), -v(16), 0]; 
            U11 = [0, 0, 0, v(20), -v(17), 0]; 
            U12 = [0, 0, 0, v(21), -v(18), 0]; 
            U13_15 = zeros(3, 6); 
            U16 = [0, 0, 0, 0, -v(10), v(7)]; 
            U17 = [0, 0, 0, 0, -v(11), v(8)]; 
            U18 = [0, 0, 0, 0, -v(12), v(9)]; 
            U19 = [0, 0, 0, v(10), 0, -v(4)]; 
            U20 = [0, 0, 0, v(11), 0, -v(5)]; 
            U21 = [0, 0, 0, v(12), 0, -v(6)]; 
            U22 = [0, 0, 0, -v(7), v(4), 0]; 
            U23 = [0, 0, 0, -v(8), v(5), 0]; 
            U24 = [0, 0, 0, -v(9), v(6), 0]; 
            dU_hat_dv = 0.5 * [U1_3; U4; U5; U6; U7; U8; U9; U10; U11; U12; U13_15; U16; U17; U18; U19; U20; U21; U22; U23; U24]; 
        end

        function S_qs = compute_S_qs(obj, n_gauss_points)
            S_qs = zeros(24, 6); 
            dx = obj.L0/2; 
            [weights, points] = obj.gauss_quadrature(n_gauss_points); 
            for i=1:length(points)
                weight = weights(i); 
                point = points(i); 
                D_gamma = obj.compute_D_gamma(point); 
                D_omega = obj.compute_D_omega(point); 
                v = D_gamma * obj.chi; 
                v2 = D_omega * obj.chi;  
                U = obj.compute_dU_bar_dv(v); 
                U2 = obj.compute_dU_hat_dv(v2); 
                S_qs = S_qs + (transpose(D_gamma) * U + transpose(D_omega) * U2) * (dx * weight); % 24 x 6 matrix. 
            end
        end


        function S_q_lambda = compute_S_q_lambda(obj)
            U13 = zeros(3, 12); % first three rows and row 13-15. 
            Z6 = zeros(1,6); 
            U4 = [obj.chi(4), 0, 0, 0, obj.chi(10), obj.chi(7), Z6]; 
            U5 = [obj.chi(5), 0, 0, 0, obj.chi(11), obj.chi(8), Z6];  
            U6 = [obj.chi(6), 0, 0, 0, obj.chi(12), obj.chi(9), Z6]; 
            U7 = [0, obj.chi(7), 0, obj.chi(10), 0, obj.chi(4), Z6]; 
            U8 = [0, obj.chi(8), 0, obj.chi(11), 0, obj.chi(5), Z6]; 
            U9 = [0, obj.chi(9), 0, obj.chi(12), 0, obj.chi(6), Z6]; 
            U10 = [0, 0, obj.chi(10), obj.chi(7), obj.chi(4), 0, Z6]; 
            U11 = [0, 0, obj.chi(11), obj.chi(8), obj.chi(5), 0, Z6]; 
            U12 = [0, 0, obj.chi(12), obj.chi(9), obj.chi(6), 0, Z6]; 
            U16 = [Z6, obj.chi(16), 0, 0, 0, obj.chi(22), obj.chi(19)]; 
            U17 = [Z6, obj.chi(17), 0, 0, 0, obj.chi(23), obj.chi(20)]; 
            U18 = [Z6, obj.chi(18), 0, 0, 0, obj.chi(24), obj.chi(21)]; 
            U19 = [Z6, 0, obj.chi(19), 0, obj.chi(22), 0, obj.chi(16)]; 
            U20 = [Z6, 0, obj.chi(20), 0, obj.chi(23), 0, obj.chi(17)];
            U21 = [Z6, 0, obj.chi(21), 0, obj.chi(24), 0, obj.chi(18)];
            U22 = [Z6, 0, 0, obj.chi(22), obj.chi(19), obj.chi(16), 0]; 
            U23 = [Z6, 0, 0, obj.chi(23), obj.chi(20), obj.chi(17), 0]; 
            U24 = [Z6, 0, 0, obj.chi(24), obj.chi(21), obj.chi(18), 0]; 
            S_q_lambda = [U13; U4; U5; U6; U7; U8; U9; U10; U11; U12; U13; U16; U17; U18; U19; U20; U21; U22; U23; U24]; 
        end             

        % returns S_q_chi matrix. 
        function S_q_chi = compute_S_q_chi(obj, n_gauss_points, C, gamma_ref, omega_ref)
            Kg = obj.compute_U2_integral_2(n_gauss_points, C, gamma_ref, omega_ref); 
            V = obj.compute_V_tot(); 
            S_q_chi = transpose(Kg) + transpose(V); 
        end

        function S_q_mu = compute_S_q_mu(obj)
            S_q_mu = -transpose(obj.compute_B(0));
        end
        
        % returns the S_q_gamma matrix.  
        function S_q_gamma = compute_S_q_gamma(obj)
            H = obj.compute_H_tot(); 
            S_q_gamma = transpose(H); 
        end

        % returns the S_q_xi matrix - constitutive eq.  
        function S_q_xi = compute_S_q_xi(obj)
            S_q_xi = zeros(24, 6); % 24x6 zero matrix. 
        end

        function S_eq = compute_S_eq(obj)
            S_eq = zeros(6, 24); % assuming dim(q) = 2 * 12. 
        end

        function S_ee = compute_S_ee(obj, C)
            % assume for now C = identity matrix. 
            S_ee = C;  
        end

        function S_es = compute_S_es(obj)
            S_es = zeros(6,6); % assuming dim(e) = dim(s) = 6. 
        end

        function S_e_lambda = compute_S_e_lambda(obj)
            S_e_lambda = zeros(6,12); 
        end

        function S_e_chi = compute_S_e_chi(obj)
            S_e_chi = zeros(6,24); 
        end

        function S_e_mu = compute_S_e_mu(obj)
            S_e_mu = eye(6); 
        end

        function S_e_gamma = compute_S_e_gamma(obj)
            S_e_gamma = zeros(6, 12); 
        end
        
        % derivative of delta s with respect to xi lagrange multipliers. 
        function S_e_xi = compute_S_e_xi(obj)
            S_e_xi = zeros(6,6); 
        end
        
        % computes the derivative of delta s with respect to q. 
        function S_sq = compute_S_sq(obj, n_gauss_points)
            [weights, points] = obj.gauss_quadrature(n_gauss_points); 
            dx = obj.L0 / 2; 
            S_sq = zeros(6, 24); 
            for i=1:length(points)
                weight = weights(i);
                point = points(i); 
                dB1_dq = obj.compute_dB1_dq(point, obj.chi); 
                dB2_dq = obj.compute_dB2_dq(point, obj.chi); 
                S_sq = S_sq + (dB1_dq + dB2_dq) * (weight * dx);  
            end
        end
        
        function S_se = compute_S_se(obj)
            S_se = zeros(6,6); 
        end

        function S_ss = compute_S_ss(obj, C_inv)
            S_ss = C_inv; % assuming C^-1 = Identity. 
        end

        function S_s_lambda = compute_S_s_lambda(obj)
            S_s_lambda = zeros(6, 12); 
        end

        function S_s_chi = compute_S_s_chi(obj, n_gauss_points)
            S_s_chi = zeros(6, 24); 
            [weights, points] = obj.gauss_quadrature(n_gauss_points); 
            dx = obj.L0 / 2; 
            for i=1:length(weights)
                weight = weights(i); 
                point = points(i); 
                B = obj.compute_B(point); 
                S_s_chi = S_s_chi + B * (weight * dx); 
            end
        end

        function S_s_mu = compute_S_s_mu(obj)
            S_s_mu = zeros(6, 6); 
        end
        
        % computes the derivative of delta s w.r.t gamma lagrange
        % multipliers. 
        function S_s_gamma = compute_S_s_gamma(obj)
            S_s_gamma = zeros(6,12); 
        end
        
        function S_s_xi = compute_S_s_xi(obj)
            S_s_xi = zeros(6, 6);
        end

        function S_lambda_q = compute_S_lambda_q(obj)            
            Z1_12 = zeros(1, 12); % 1x12 zero matrix. 
            Z13 = zeros(1,3); % 1x3 zero matrix. 
            chi1 = obj.chi(1:3); 
            chi2 = obj.chi(4:6); 
            chi3 = obj.chi(7:9); 
            chi4 = obj.chi(10:12); 
            chi5 = obj.chi(13:15); 
            chi6 = obj.chi(16:18); 
            chi7 = obj.chi(19:21); 
            chi8 = obj.chi(22:24); 
            dH_Lambda_dq_1 = [Z13, transpose(chi2), Z13, Z13, Z1_12]; 
            dH_Lambda_dq_2 = [Z13, Z13, transpose(chi3), Z13, Z1_12]; 
            dH_Lambda_dq_3 = [Z13, Z13, Z13, transpose(chi4), Z1_12]; 
            dH_Lambda_dq_4 = [Z13, Z13, transpose(chi4), transpose(chi3), Z1_12]; 
            dH_Lambda_dq_5 = [Z13, transpose(chi4), Z13, transpose(chi2), Z1_12]; 
            dH_Lambda_dq_6 = [Z13, transpose(chi3), transpose(chi2), Z13, Z1_12]; 
            dH_Lambda_dq_7 = [Z1_12, Z13, transpose(chi5), Z13, Z13]; 
            dH_Lambda_dq_8 = [Z1_12, Z13, Z13, transpose(chi7), Z13]; 
            dH_Lambda_dq_9 = [Z1_12, Z13, Z13, Z13, transpose(chi8)]; 
            dH_Lambda_dq_10 = [Z1_12, Z13, Z13, transpose(chi8), transpose(chi7)];
            dH_Lambda_dq_11 = [Z1_12, Z13, transpose(chi8), Z13, transpose(chi6)];
            dH_Lambda_dq_12 = [Z1_12, Z13, transpose(chi7), transpose(chi6), Z13]; 
            S_lambda_q = [dH_Lambda_dq_1; dH_Lambda_dq_2; dH_Lambda_dq_3; dH_Lambda_dq_4; dH_Lambda_dq_5; dH_Lambda_dq_6; dH_Lambda_dq_7;dH_Lambda_dq_8; dH_Lambda_dq_9; dH_Lambda_dq_10; dH_Lambda_dq_11; dH_Lambda_dq_12];
        end

        function S_lambda_e = compute_S_lambda_e(obj)
            S_lambda_e = zeros(12, 6); % 12x6 zero matrix. 
        end

        function S_lambda_s = compute_S_lambda_s(obj)
            S_lambda_s = zeros(12, 6); % 12x6 zero matrix. 
        end

        function S_lambda_lambda = compute_S_lambda_lambda(obj)
            S_lambda_lambda = zeros(12,12); % 12x12 zero matrix. 
        end
        
        function S_lambda_chi = compute_S_lambda_chi(obj)
            S_lambda_chi = obj.compute_H_tot(); % 12x24 matrix. 
        end

        function S_lambda_gamma = compute_S_lambda_gamma(obj)
            S_lambda_gamma = zeros(12, 12); % 12x12 zero matrix. 
        end

        function S_lambda_mu = compute_S_lambda_mu(obj)
            S_lambda_mu = zeros(12,6); % 12x6 zero matrix. 
        end
        
        function S_lambda_xi = compute_S_lambda_xi(obj)
            S_lambda_xi = zeros(12,6); % 12x6 zero matrix.          
        end

        function S_chi_q = compute_S_chi_q(obj, n_gauss_points, C)
            U2 = obj.compute_U2_integral_2(n_gauss_points, C, obj.gamma_ref, obj.omega_ref);
            V = obj.compute_V_tot(); 
            S_chi_q = U2 + V; 
        end

        function S_chi_e = compute_S_chi_e(obj)
            S_chi_e = zeros(24, 6); 
        end

        function S_chi_s = compute_S_chi_s(obj, n_gauss_points) 
            [weights, points] = obj.gauss_quadrature(n_gauss_points); 
            dx = obj.L0/2; 
            S_chi_s = zeros(24, 6); 
            for i=1:length(weights)
                weight = weights(i); 
                point = points(i); 
                B = obj.compute_B(point); 
                S_chi_s = S_chi_s + transpose(B) * (weight * dx); 
            end
        end

        function S_chi_lambda = compute_S_chi_lambda(obj)
            H = obj.compute_H_tot(); 
            S_chi_lambda = transpose(H); 
        end

        function S_chi_chi = compute_S_chi_chi(obj)
            S_chi_chi = zeros(24, 24); 
        end

        function S_chi_mu = compute_S_chi_mu(obj)
            S_chi_mu = zeros(24, 6); 
        end

        function S_chi_gamma = compute_S_chi_gamma(obj)
            S_chi_gamma = zeros(24, 12); 
        end

        function S_chi_xi = compute_S_chi_xi(obj)
            S_chi_xi = zeros(24, 6); 
        end

        function S_mu_q = compute_S_mu_q(obj)
            S_mu_q = -obj.compute_B(0); % 6x24
        end

        function S_mu_e = compute_S_mu_e(obj)
            S_mu_e = eye(6); 
        end

        function S_mu_s = compute_S_mu_s(obj)
            S_mu_s = zeros(6, 6); 
        end

        function S_mu_lambda = compute_S_mu_lambda(obj)
            S_mu_lambda = zeros(6, 12);
        end

        function S_mu_chi = compute_S_mu_chi(obj)
            S_mu_chi = zeros(6, 24); 
        end
        
        function S_mu_mu = compute_S_mu_mu(obj)
            S_mu_mu = zeros(6,6);
        end

        function S_mu_gamma = compute_S_mu_gamma(obj)
            S_mu_gamma = zeros(6, 12); 
        end

        function S_mu_xi = compute_S_mu_xi(obj)
            S_mu_xi = zeros(6, 6);
        end
        
        function S_gamma_q = compute_S_gamma_q(obj)
            S_gamma_q = obj.compute_H_tot(); % 12 x 24 matrix. 
        end

        function S_gamma_e = compute_S_gamma_e(obj)
            S_gamma_e = zeros(12, 6); % 12 x 24 zero matrix.  
        end

        function S_gamma_s = compute_S_gamma_s(obj)
            S_gamma_s = zeros(12, 6); % 12 x 6 zero matrix.  
        end

        function S_gamma_lambda = compute_S_gamma_lambda(obj)
            S_gamma_lambda = zeros(12, 12); % 12 x 12 zero matrix.
        end

        function S_gamma_chi = compute_S_gamma_chi(obj)
            S_gamma_chi = zeros(12, 24); % 12 x 24 zero matrix. 
        end

        function S_gamma_mu = compute_S_gamma_mu(obj)
            S_gamma_mu = zeros(12, 6); % 12 x 6 zero matrix. 
        end

        function S_gamma_gamma = compute_S_gamma_gamma(obj)
            S_gamma_gamma = zeros(12, 12); % 12x12 zero matrix.
        end

        function S_gamma_xi = compute_S_gamma_xi(obj)
            S_gamma_xi = zeros(12, 6); % 12x6 zero matrix. 
        end

        function S_xi_q = compute_S_xi_q(obj)
            S_xi_q = zeros(6, 24); % 6 x 24 zero matrix. 
        end

        function S_xi_e = compute_S_xi_e(obj)
            S_xi_e = zeros(6,6); % 6x6 zero matrix. 
        end

        function S_xi_s = compute_S_xi_s(obj)
            S_xi_s = zeros(6,6); % 6x6 zero matrix. 
        end

        function S_xi_lambda = compute_S_xi_lambda(obj)
            S_xi_lambda = zeros(6,12); % 6x12 zero matrix. 
        end

        function S_xi_chi = compute_S_xi_chi(obj)
            S_xi_chi = zeros(6, 24); % 6x24 zero matrix. 
        end

        function S_xi_mu = compute_S_xi_mu(obj)
            S_xi_mu = zeros(6, 6); % 6x6 zero matrix. 
        end

        function S_xi_gamma = compute_S_xi_gamma(obj)
            S_xi_gamma = zeros(6,12); % 6x12 zero matrix.  
        end

        function S_xi_xi = compute_S_xi_xi(obj)
            S_xi_xi = zeros(6, 6); % 6x6 zero matrix. 
        end

        % assembling the submatrix S_fix for the KKT matrix. 
        function S_fix = compute_KKT_S_fix(obj, C_material, C, n_gauss_points)
            % compute S_qq. 
            S_qq = obj.compute_S_qq(); 
            % compute S_qe. 
            S_qe = zeros(24, 6); 
            % compute S_qs. 
            S_qs = obj.compute_S_qs(n_gauss_points); 
            % compute S_q_lambda. 
            S_q_lambda = obj.compute_S_q_lambda(); 
            % compute S_q_chi. 
            S_q_chi = obj.compute_S_q_chi(n_gauss_points, C_material, obj.gamma_ref, obj.omega_ref); 
            % compute S_q_mu. 
            S_q_mu = obj.compute_S_q_mu(); 
            % compute S_gamma. 
            S_q_gamma = obj.compute_S_q_gamma();           
            S1 = [S_qq, S_qe, S_qs, S_q_lambda, S_q_chi, S_q_mu, S_q_gamma]; 
            S2 = [obj.compute_S_eq(), obj.compute_S_ee(C), obj.compute_S_es(), obj.compute_S_e_lambda(), obj.compute_S_e_chi(), obj.compute_S_e_mu(), obj.compute_S_e_gamma()]; 
            S3 = [obj.compute_S_sq(n_gauss_points), obj.compute_S_se(), obj.compute_S_ss(C), obj.compute_S_s_lambda(), obj.compute_S_s_chi(n_gauss_points), obj.compute_S_s_mu(), obj.compute_S_s_gamma()]; 
            S4 = [obj.compute_S_lambda_q(), obj.compute_S_lambda_e(), obj.compute_S_lambda_s(), obj.compute_S_lambda_lambda(), obj.compute_S_lambda_chi(), obj.compute_S_lambda_mu(), obj.compute_S_lambda_gamma()]; 
            S5 = [obj.compute_S_chi_q(n_gauss_points, C_material), obj.compute_S_chi_e(), obj.compute_S_chi_s(n_gauss_points), obj.compute_S_chi_lambda(), obj.compute_S_chi_chi(), obj.compute_S_chi_mu(), obj.compute_S_chi_gamma()];
            S6 = [obj.compute_S_mu_q(), obj.compute_S_mu_e(), obj.compute_S_mu_s(), obj.compute_S_mu_lambda(), obj.compute_S_mu_chi(), obj.compute_S_mu_mu(), obj.compute_S_mu_gamma()]; 
            S7 = [obj.compute_S_gamma_q(), obj.compute_S_gamma_e(), obj.compute_S_gamma_s(), obj.compute_S_gamma_lambda(), obj.compute_S_gamma_chi(), obj.compute_S_gamma_mu(), obj.compute_S_gamma_gamma()]; 
            S_fix = [S1; S2; S3; S4; S5; S6; S7];
        end

        % assembles the S matrix for KKT. 
        function S = compute_KKT_S(obj, C_material, C, n_gauss_points, e_tilde, s_tilde)
            % S_e_tilde_e_tilde. 
            S_e_tilde_e_tilde = obj.compute_S_e_tilde_e_tilde(C); 
            % S_e_tilde_s_tilde. 
            S_e_tilde_s_tilde = obj.compute_S_e_tilde_s_tilde(C); 
            S_e_tilde_x = obj.compute_S_e_tilde_x(C); 
            S_e_tilde_xi = obj.compute_S_e_tilde_xi(C, e_tilde, s_tilde); 
            % row 1. 
            S1 = [S_e_tilde_e_tilde, S_e_tilde_s_tilde, S_e_tilde_x, S_e_tilde_xi]; 
            S_s_tilde_e_tilde = obj.compute_S_s_tilde_e_tilde(); % 6x6
            S_s_tilde_s_tilde = obj.compute_S_s_tilde_s_tilde(C); % 6x6
            S_s_tilde_x = obj.compute_S_s_tilde_x(C); % 6x90 
            S_s_tilde_xi = obj.compute_S_s_tilde_xi(C, e_tilde, s_tilde); % 6x6
            % row 2. 
            S2 = [S_s_tilde_e_tilde, S_s_tilde_s_tilde, S_s_tilde_x, S_s_tilde_xi]; 
            S_q_e_tilde = zeros(24, 6); 
            S_q_s_tilde = zeros(24, 6); 
            S_qq = obj.compute_S_qq(); % 24x24
            S_qe = obj.compute_S_qe(); % 24x6
            S_qs = obj.compute_S_qs(n_gauss_points); % 24x6 
            S_q_lambda = obj.compute_S_q_lambda(); % 21 x 12 -> check
            S_q_chi = obj.compute_S_q_chi(n_gauss_points, C_material, obj.gamma_ref, obj.omega_ref); 
            S_q_mu = obj.compute_S_q_mu(); 
            S_q_gamma = obj.compute_S_q_gamma(); 
            S_q_xi = obj.compute_S_q_xi();
            % row 3. 
            S3 = [S_q_e_tilde, S_q_s_tilde, S_qq, S_qe, S_qs, S_q_lambda, S_q_chi, S_q_mu, S_q_gamma, S_q_xi]; 
            S_e_e_tilde = -transpose(C); % 6x6
            S_e_s_tilde = zeros(6,6); % 6x6
            S_eq = obj.compute_S_eq(); % 6x24
            S_ee = obj.compute_S_ee(C); % 6x6
            S_es = obj.compute_S_es(); % 6x6
            S_e_lambda = obj.compute_S_e_lambda(); % 6x12
            S_e_chi = obj.compute_S_e_chi(); % 6x24
            S_e_mu = obj.compute_S_e_mu(); % 6x6
            S_e_gamma = obj.compute_S_e_gamma(); 
            S_e_xi = obj.compute_S_e_xi(); 
            % row 4.
            S4 = [S_e_e_tilde, S_e_s_tilde, S_eq, S_ee, S_es, S_e_lambda, S_e_chi, S_e_mu, S_e_gamma, S_e_xi]; 
            S_s_e_tilde = zeros(6,6); % 6x6
            S_s_s_tilde = -transpose(C); % 6x6
            S_sq = obj.compute_S_sq(n_gauss_points);
            S_se = obj.compute_S_se(); 
            S_ss = obj.compute_S_ss(C); 
            S_s_lambda = obj.compute_S_s_lambda(); % 6x12
            S_s_chi = obj.compute_S_s_chi(n_gauss_points); % 6x24
            S_s_mu = obj.compute_S_s_mu(); % 6x6
            S_s_gamma = obj.compute_S_s_gamma(); 
            S_s_xi = obj.compute_S_s_xi(); 
            % row 5. 
            S5 = [S_s_e_tilde, S_s_s_tilde, S_sq, S_se, S_ss, S_s_lambda, S_s_chi, S_s_mu, S_s_gamma, S_s_xi]; 
            S_lambda_e_tilde = zeros(12, 6); % 12x6
            S_lambda_s_tilde = zeros(12, 6); % 12x6
            S_lambda_q = obj.compute_S_lambda_q(); 
            S_lambda_e = obj.compute_S_lambda_e(); 
            S_lambda_s = obj.compute_S_lambda_s(); 
            S_lambda_lambda = obj.compute_S_lambda_lambda(); 
            S_lambda_chi = obj.compute_S_lambda_chi(); 
            S_lambda_mu = obj.compute_S_lambda_mu();
            S_lambda_gamma = obj.compute_S_lambda_gamma(); 
            S_lambda_xi = obj.compute_S_lambda_xi(); 
            % row 6. 
            S6 = [S_lambda_e_tilde, S_lambda_s_tilde, S_lambda_q, S_lambda_e, S_lambda_s, S_lambda_lambda, S_lambda_chi, S_lambda_mu, S_lambda_gamma, S_lambda_xi]; 
            S_chi_e_tilde = zeros(24, 6); 
            S_chi_s_tilde = zeros(24, 6); 
            S_chi_q = obj.compute_S_chi_q(n_gauss_points, C_material); 
            S_chi_e = obj.compute_S_chi_e(); 
            S_chi_s = obj.compute_S_chi_s(n_gauss_points); 
            S_chi_lambda = obj.compute_S_chi_lambda(); 
            S_chi_chi = obj.compute_S_chi_chi(); 
            S_chi_mu = obj.compute_S_chi_mu(); % 24x6
            S_chi_gamma = obj.compute_S_chi_gamma(); % 24x12 
            S_chi_xi = obj.compute_S_chi_xi(); % 24x6
            % row 7. 
            S7 = [S_chi_e_tilde, S_chi_s_tilde, S_chi_q, S_chi_e, S_chi_s, S_chi_lambda, S_chi_chi, S_chi_mu, S_chi_gamma, S_chi_xi];            
            S_mu_e_tilde = zeros(6, 6); 
            S_mu_s_tilde = zeros(6, 6); 
            S_mu_q = obj.compute_S_mu_q(); 
            S_mu_e = obj.compute_S_mu_e(); 
            S_mu_s = obj.compute_S_mu_s();
            S_mu_lambda = obj.compute_S_mu_lambda(); 
            S_mu_chi = obj.compute_S_mu_chi(); 
            S_mu_mu = zeros(6, 6); 
            S_mu_gamma = obj.compute_S_mu_gamma(); 
            S_mu_xi = obj.compute_S_mu_xi(); 
            % row 8. 
            S8 = [S_mu_e_tilde, S_mu_s_tilde, S_mu_q, S_mu_e, S_mu_s, S_mu_lambda, S_mu_chi, S_mu_mu, S_mu_gamma, S_mu_xi]; 
            S_gamma_e_tilde = zeros(12,6); 
            S_gamma_s_tilde = zeros(12,6); 
            S_gamma_q = obj.compute_S_gamma_q(); 
            S_gamma_e = obj.compute_S_gamma_e(); 
            S_gamma_s = obj.compute_S_gamma_s(); % 12x6
            S_gamma_lambda = obj.compute_S_gamma_lambda(); % 12x12
            S_gamma_chi = obj.compute_S_gamma_chi(); % 12x24
            S_gamma_mu = obj.compute_S_gamma_mu(); % 12x6
            S_gamma_gamma = obj.compute_S_gamma_gamma(); % 12x12
            S_gamma_xi = obj.compute_S_gamma_xi(); % 12x6
            % row 9. 
            S9 = [S_gamma_e_tilde, S_gamma_s_tilde, S_gamma_q, S_gamma_e, S_gamma_s, S_gamma_lambda, S_gamma_chi, S_gamma_mu, S_gamma_gamma, S_gamma_xi]; 
            S_xi_e_tilde = obj.compute_dg_de_tilde(C, e_tilde, s_tilde); % zeros(6,6);     
            S_xi_s_tilde = obj.compute_dg_ds_tilde(C, e_tilde, s_tilde); 
            S_xi_q = obj.compute_S_xi_q(); % 6x24
            S_xi_e = obj.compute_S_xi_e(); % 6x6
            S_xi_s = obj.compute_S_xi_s(); % 6x6
            S_xi_lambda = obj.compute_S_xi_lambda(); % 6x12
            S_xi_chi = obj.compute_S_xi_chi(); % 6x24
            S_xi_mu = obj.compute_S_xi_mu(); % 6x6
            S_xi_gamma = obj.compute_S_xi_gamma(); % 6x12
            S_xi_xi = obj.compute_S_xi_xi(); % 6x6
            % row 10. 
            S10 = [S_xi_e_tilde, S_xi_s_tilde, S_xi_q, S_xi_e, S_xi_s, S_xi_lambda, S_xi_chi, S_xi_mu, S_xi_gamma, S_xi_xi]; 
            % S = [S1; S2; S3; S4; S5; S6; S7; S8; S9; S10]; 
            S_fix = obj.compute_KKT_S_fix(C_material, C, n_gauss_points); 
            row_1 = [S_e_tilde_e_tilde, S_e_tilde_s_tilde, S_e_tilde_x, S_e_tilde_xi]; 
            row_2 = [transpose(S_e_tilde_s_tilde), S_s_tilde_s_tilde, S_s_tilde_x, S_s_tilde_xi]; 
            row_3 = [transpose(S_e_tilde_x), transpose(S_s_tilde_x), S_fix, zeros(90, 6)]; 
            row_4 = [transpose(S_e_tilde_xi), transpose(S_s_tilde_xi), zeros(6, 90), zeros(6, 6)];
            S = [row_1; row_2; row_3; row_4]; 
        end

        % first order KKT conditions.
        % variation of the Lagrangian w.r.t q. 
        function dL_dq = compute_dL_dq(obj, n_gauss_points, C)
            % computes the derivative of integral B^T * s w.r.t q.
            dBs_dq = obj.compute_U2_integral_2(n_gauss_points, C, obj.gamma_ref, obj.omega_ref); 
            % computes the derivative of H^T * lambda w.r.t q. 
            V = obj.compute_V_tot(); 
            % compute B(q). 
            B = obj.compute_B(0); % s = 0 interpolation value. 
            % compute H(q). 
            H = obj.compute_H_tot();
            % compute dL_dq. 
            dL_dq = transpose(obj.chi) * (dBs_dq + V) - transpose(obj.mu) * B + transpose(obj.gamma) * H; 
        end

        % variation of the Lagrangian w.r.t e. 
        % C - coefficient matrix.
        function dL_de = compute_dL_de(obj, e_tilde, C)
            delta_e = obj.e0 - e_tilde; 
            dL_de = transpose(delta_e) * C + transpose(obj.mu);  
        end

        % variation of the Lagrangian w.r.t s. 
        function dL_ds = compute_dL_ds(obj, s_tilde, C, n_gauss_points)
            [weights, points] = obj.gauss_quadrature(n_gauss_points); 
            % compute the integral of B^T(q).
            Bt_q = 0;
            dx = obj.L0/2; 
            for i = 1:length(weights)
                weight = weights(i); 
                point = points(i); 
                % compute B. 
                B = obj.compute_B(point); 
                Bt_q = Bt_q + transpose(B) * (weight * dx); 
            end
            % compute dL_ds. 
            delta_s = obj.s0 - s_tilde; 
            dL_ds = transpose(delta_s) * C + transpose(obj.chi) * Bt_q; % 1x6 vector. 
        end

        % variation of the Lagrangian w.r.t lambda. 
        function dL_dLamba = compute_dL_dLambda(obj)
            H = obj.compute_H_tot(); 
            dL_dLamba = transpose(obj.chi) * transpose(H); % 1x12 vector. 
        end

        % variation of the Lagrangian w.r.t chi. 
        function dL_dChi = compute_dL_dChi(obj, n_gauss_points, C, f_ext)
            f_int = obj.compute_f_int(n_gauss_points, C, obj.gamma_ref, obj.omega_ref); % internal forces. 
            H = obj.compute_H_tot();
            f_H = transpose(H) * [obj.lambda_At; obj.lambda_Bt]; % constraint forces. 
            f_balance = f_int + f_H - f_ext; % balance eq. 
            dL_dChi = transpose(f_balance); 
        end

        % variation of the Lagrangian w.r.t muh. 
        function dL_dmuh = compute_dL_dmuh(obj)
            e_q = obj.compute_e(0, obj.gamma_ref, obj.omega_ref); 
            delta_e = obj.e0 - e_q; 
            dL_dmuh = transpose(delta_e); 
        end

        % variation of the Lagrangian w.r.t gamma coeffs. 
        function dL_dgamma = compute_dL_dgamma(obj)
            h = obj.compute_h_q(); % h(q) vector - restrictions condition. 
            dL_dgamma = transpose(h); 
        end
        
        % variation of the Lagrangian w.r.t xi coeffs. 
        function dL_dxi = compute_dL_dxi(obj, C, e_tilde, s_tilde)
            % the constitutive eq. 
            g = obj.compute_constitutive_g(C, e_tilde, s_tilde); 
            dL_dxi = transpose(g); 
        end

        % variation of the Lagrangian w.r.t e_tilde. 
        function dL_de_tilde = compute_dL_de_tilde(obj, C_material, C, e_tilde, s_tilde)
            % compute the derivative of the cost function. 
            delta_e = obj.e0 - e_tilde; 
            % compute the derivative of the constitutive eq w.r.t e_tilde. 
            dg_de_tilde = obj.compute_dg_de_tilde(C_material, e_tilde, s_tilde); 
            dL_de_tilde = -transpose(delta_e) * C + transpose(obj.xi) * dg_de_tilde; % 1x6 vector.
        end

        % variation of the Lagrangian w.r.t s_tilde. 
        function dL_ds_tilde = compute_dL_ds_tilde(obj, C_material, C, e_tilde, s_tilde)
            delta_s = obj.s0 - s_tilde;
            % compute the derivative of the constitutive eq w.r.t s_tilde. 
            dg_ds_tilde = obj.compute_dg_ds_tilde(C_material, e_tilde, s_tilde); 
            dL_ds_tilde = -transpose(delta_s) * C  + transpose(obj.xi) * dg_ds_tilde; % 1x6 vector. 
           
        end

        % computes the derivative vector of the lagrangian.
        % C_material - related to the constitutive eq. 
        % C - for computing cost function. 
        function dL = compute_dL(obj, C_material, C, e_tilde, s_tilde, n_gauss_points, f_ext)
            % compute dL_de_tilde. 
            dL_de_tilde = obj.compute_dL_de_tilde(C_material, C, e_tilde, s_tilde);
            % compute dL_ds_tilde. 
            dL_ds_tilde = obj.compute_dL_ds_tilde(C_material, C, e_tilde, s_tilde); 
            % compute dL_dq. 
            dL_dq = obj.compute_dL_dq(n_gauss_points, C_material); % 1x24
            % compute dL_de. 
            dL_de = obj.compute_dL_de(e_tilde, C); % 1x6
            % compute dL_ds. 
            dL_ds = obj.compute_dL_ds(s_tilde, C, n_gauss_points); 
            % compute dL_dLambda. 
            dL_dLambda = obj.compute_dL_dLambda(); 
            % compute dL_dChi. 
            dL_dChi = obj.compute_dL_dChi(n_gauss_points, C_material, f_ext); 
            % compute dL_dmu. 
            dL_dmu = obj.compute_dL_dmuh(); 
            % compute dL_dgamma. 
            dL_dgamma = obj.compute_dL_dgamma(); 
            % compute dL_dxi. 
            dL_dxi = obj.compute_dL_dxi(C_material, e_tilde, s_tilde); 
            % dL vector. 
            dL = [dL_de_tilde, dL_ds_tilde, dL_dq, dL_de, dL_ds, dL_dLambda, dL_dChi, dL_dmu, dL_dgamma, dL_dxi]; 
            dL = transpose(dL); % 108 x 1 vector. 
        end

        % for computing the cost function of the Lagrangian.
        function J_cost = compute_J_cost(obj, e_tilde, s_tilde)
            delta_e = obj.e0 - e_tilde; 
            delta_s = obj.s0 - s_tilde; 
            % assume C = identity matrix. 
            J_cost = dot(delta_e, delta_e) * (1/2) + dot(delta_s, delta_s) * (1/2); 
        end
        
        % for computing the compatibility eq. mu * (e - e(q))
        function compatibilty = compute_compatibility_constraint(obj)
            e_q = obj.compute_e(); 
            delta_e = obj.e0 - e_q; 
            compatibilty = transpose(obj.mu) * delta_e; 
        end
        
        % for computing the balance constraint eq. 
        function balance = compute_balance_constraint(obj, n_gauss_points, C, f_ext)           
            % compute the internal force. 
            f_int = obj.compute_f_int(n_gauss_points, C, obj.gamma_ref, obj.omega_ref); 
            % compute the constraint from H(q). 
            f_H = [obj.compute_f_H1(); obj.compute_f_H2()];
            % compute the balance. 
            balance = transpose(obj.chi) * (f_int + f_H - f_ext); 
        end

        % for computing the gamma * h(q) constraint. 
        function restriction = compute_restriction_constraint(obj)
            h = obj.compute_h_q();
            restriction = transpose(obj.gamma) * h; 
        end
        
        % for computing the constitutive equation - linear. 
        function g = compute_constitutive_g(obj, C, e_tilde, s_tilde)
            g = s_tilde - C * e_tilde; 
        end 

        % computes the derivative of the constitutive eq w.r.t e_tilde. 
        function dg_de_tilde = compute_dg_de_tilde(obj, C, e_tilde, s_tilde)
            dg_de_tilde = -C; 
        end
        
        % computes the derivative of the constitutive eq w.r.t s_tilde. 
        function dg_ds_tilde = compute_dg_ds_tilde(obj, C, e_tilde, s_tilde)
            [n,m] = size(s_tilde); 
            s = max(n,m);
            dg_ds_tilde = eye(s); 
        end

        % for computing the xi * g(e_tilde, s_tilde) constraint -
        % constitutive eq. 
        function constitutive = compute_constitutive_constraint(obj, C, e_tilde, s_tilde)
           g = obj.compute_constitutive_g(C, e_tilde, s_tilde); 
           constitutive = transpose(obj.xi) * g; 
        end

        % computes the lagrangian for a given set of strain / stress
        % variables. 
        function L_fix = compute_L_fix(obj, n_gauss_points, C, f_ext, e_tilde, s_tilde)
            % compute the cost function. 
            J_cost = obj.compute_J_cost(e_tilde, s_tilde); 
            % compute the compatibility eq. 
            compatibility_eq = obj.compute_compatibility_constraint(); 
            % compute the balance eq. 
            balance_eq = obj.compute_balance_constraint(n_gauss_points, C, f_ext); 
            % restriction eq. 
            restriction_eq = obj.compute_restriction_constraint(); 
            % constitutive eq. 
            constitutive_eq = obj.compute_constitutive_constraint(C, e_tilde, s_tilde); 
            % add all the terms. 
            L_fix = J_cost + compatibility_eq + balance_eq + restriction_eq + constitutive_eq; 
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
        
        % updates the chi Lagrange coefficients. 
        % delta_chi - 24x1 vector.
        function update_chi(obj, delta_chi)
            obj.chi = obj.chi + delta_chi; 
        end
        
        % updates the mu Lagrange coefficients. 
        function update_mu(obj, delta_mu)
            obj.mu = obj.mu + delta_mu; 
        end
        
        % updates the gamma Lagrange coefficients. 
        function update_gamma(obj, delta_gamma)
            obj.gamma = obj.gamma + delta_gamma;
        end
        
        % updates the xi Lagrange coefficients.     
        function update_xi(obj, delta_xi)
            obj.xi = obj.xi + delta_xi;
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
            



