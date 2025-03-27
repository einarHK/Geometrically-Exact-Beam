
clear;
clc; 
% comparison - single beam element with beam struct. 
%% testing
L = 1;
n_elems = 1; 
x1 = [0, 0, 0]'; 
x2 = [L, 0, 0]';
d3 = [1, 0, 0]'; % axial director. 
d2 = [0, 1, 0]'; % shear directors 2, 1. 
d1 = [0, 0, 1]'; 
dof = 12;  
n_constraints = 6;
lambda = zeros(n_constraints,1);
s = 0; 

fixed_dof = 1:12; % set first node fixed, for all kinematic variables.
fixed_const = 2*dof+1:2*dof+n_constraints; % fixed constraints at first node. 
all_dof = 1:2*dof + 2*n_constraints; % all dofs + constraints. 

free_dof = setdiff(all_dof, [fixed_dof, fixed_const]);
C = eye(6);

gamma_ref = [0;0; 1];
% gamma_ref = [0;0; 0];
omega_ref = [0;0;0];

beam = Beam_elem(x1, x2, L, d3, d2, d1, d3, d2, d1, lambda, lambda, dof, n_constraints, fixed_dof, gamma_ref, omega_ref); 
% f_int = beam.compute_f_int(1, C, gamma_ref, omega_ref); 

lambdas = repmat(lambda, 1, n_elems+1); 
beam_dofs = repmat(dof, 1, n_elems +1);
beam_constraints = repmat(n_constraints, 1, n_elems + 1); 
beam_fixed_dofs = zeros(1, n_elems + 1);
beam_fixed_dofs(1) = 12;

% calculate beam coordinates and beam lengths. 
[beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
beam_directors = repmat(vertcat(d1, d2, d3), 1, n_elems + 1);

% create beam elements and beam class object. 
beam_elems = create_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
beam2 = Beam(n_elems, beam_elems, dof, n_constraints); 
%% Single beam element - vertical load. 
Fx = 0; 
Fy = 0; 
Fz = 0.4; 

% number of load steps. 
load_steps = 10; 

% force increments. 
dFx = Fx / load_steps;
dFy = Fy / load_steps;
dFz = Fz / load_steps;

% 
n_gauss_points = 1; 

% 
max_iter = 50; 

% 
f_ext = zeros(2 * beam.dof, 1); 

% 
TOL = 1e-8;

delta_u_lambdas_beam1 = []; 

counter = 0; 

% main loop. 
for i=1:load_steps
    iter = 0; 
    F_x = dFx * i; 
    F_y = dFy * i; 
    F_z = dFz * i; 
    f_ext(beam.dof + 1) = F_x; 
    f_ext(beam.dof + 2) = F_y;
    f_ext(beam.dof + 3) = F_z; 
    while (iter < max_iter)
        iter = iter + 1; 
        counter = counter + 1; 
        % initialize solution vector. 
        delta_u_lambda = zeros(2 * beam.dof + 2 * beam.n_constraints, 1);
        
        % compute S matrix. 
        S = beam.compute_S_mat(n_gauss_points, C, gamma_ref, omega_ref); 
        S = S(free_dof, free_dof); 
        
        % compute the force terms. 
        f_int = beam.compute_f_int(n_gauss_points, C, gamma_ref, omega_ref); 
        f_H1 = beam.compute_f_H1(); % force coming from node 1 constraints. 
        f_H2 = beam.compute_f_H2(); % forces coming from node 2 constraints. 
        f_H = [f_H1; f_H2]; 
        residual = f_int + f_H - f_ext; % residual force vector. 
        
        h_q = beam.compute_h_q(); % constraint vector. 
        rhs = [-residual; -h_q]; % right hand side vector. 
        rhs = rhs(free_dof); % remove fixed indices. 

        % solve linearized system of equations. 
        delta_u_lambda(free_dof) = S \ rhs; 

        delta_u_lambdas_beam1 = [delta_u_lambdas_beam1, delta_u_lambda]; 

        % extract the beam value increments. 
        dx1 = delta_u_lambda(1:3); 
        delta_d1_A = delta_u_lambda(4:6); 
        delta_d2_A = delta_u_lambda(7:9); 
        delta_d3_A = delta_u_lambda(10:12); 
        dx2 = delta_u_lambda(13:15); 
        delta_d1_B = delta_u_lambda(16:18); 
        delta_d2_B = delta_u_lambda(19:21); 
        delta_d3_B = delta_u_lambda(22:24); 
        delta_lambda_A = delta_u_lambda(25:24 + n_constraints, 1); 
        delta_lambda_B = delta_u_lambda(25 + n_constraints: 24 + 2*n_constraints, 1); 

        % updating the beam parameters. 
        beam.update_params(dx1, dx2, delta_d1_A, delta_d2_A, delta_d3_A, delta_d1_B, delta_d2_B, delta_d3_B);
        beam.update_lambda(delta_lambda_A, delta_lambda_B); 

        % compute the norm of the rhs vector. 
        norm_error = norm(rhs); 

        % if the norm less than tolerance, break out of loop. 
        if (norm_error < TOL)
            break;
        end

    end
end 

% display end node position and directors. 
beam.display_node_pos(); 
beam.display_node_directors();
%% Beam struct - vertical force. 
Fx = 0; 
Fy = 0; 
Fz = 0.4; 

% number of load steps. 
load_steps = 10; 

% force increments. 
dFx = Fx / load_steps;
dFy = Fy / load_steps;
dFz = Fz / load_steps;

delta_u_lambdas_beam2 = []; 

% 
n_gauss_points = 1; 

% 
max_iter = 50; 

% 
TOL = 1e-8;

free_dof = beam2.compute_free_dof();
f_ext = zeros(beam2.dof_per_node * beam2.n_nodes, 1); 

% main loop. 
for i=1:load_steps
    iter = 0; 
    F_x = dFx * i; 
    F_y = dFy * i; 
    F_z = dFz * i; 
    f_ext(beam.dof * (n_elems) + 1) = F_x; 
    f_ext(beam.dof * (n_elems) + 2) = F_y;
    f_ext(beam.dof * (n_elems) + 3) = F_z; 
    while (iter < max_iter)
        iter = iter + 1; 
        % initialize solution vector. 
        delta_u_lambda = zeros(beam2.n_nodes * (beam2.dof_per_node + beam2.constraint_per_node), 1);
        
        % compute S matrix. 
        S = beam2.compute_S_mat(n_gauss_points, C); 
        S = S(free_dof, free_dof); 
        
        % compute the force terms. 
        f_int = beam2.compute_f_int(n_gauss_points, C);
        f_H = beam2.compute_f_H();
        residual = f_int + f_H - f_ext; % residual force vector. 
        
        h_q = beam2.compute_h_q(); % constraint vector. 
        rhs = [-residual; -h_q]; % right hand side vector. 
        rhs = rhs(free_dof); % remove fixed indices. 

        % solve linearized system of equations. 
        delta_u_lambda(free_dof) = S \ rhs; 
        
        delta_u_lambdas_beam2 = [delta_u_lambdas_beam2, delta_u_lambda]; 

        % update the beam element parameters. 
        beam2.update_beam_params(delta_u_lambda ); 

        % compute the norm of the rhs vector. 
        norm_error = norm(rhs); 

        % if the norm less than tolerance, break out of loop. 
        if (norm_error < TOL)
            break;
        end

    end
end

% display end node position and directors. 
beam2.display_end_node_directors(); 
beam2.display_end_node_pos();
 
beam2.plot(7, 7, 7, 0, "");

%% 

