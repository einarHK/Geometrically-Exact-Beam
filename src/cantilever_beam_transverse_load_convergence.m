% cantilever beam - transverse loads, convergence plot as function of beam 
% elements. 
clear; 
clc; 

%% CONSTANTS. 
n_elems = [2, 4, 8, 16]; % number of beam elements for each case. 
L0 = 102.75; % initial length of the beam. 

x1 = [0, 0, 0]'; % start position of the beam.
x2 = [L0, 0, 0]'; % end position of the beam. 
d3 = [1, 0, 0]'; % axial director. 
d2 = [0, 1, 0]'; % shear directors 2, 1. 
d1 = [0, 0, 1]'; 
dof = 12; % kinematic degrees of freedom per node. 
fixed_dof = 12; % fix the kinematic dof at the first node. 
n_constraints = 6; % number of constraints per node. 
lambda = zeros(n_constraints,1); % lambda multipliers per node - initialized as zeros. 

% discretize the domain into beam coordinates, beam lengths and directors. 
all_beams = []; 
for i=1:length(n_elems)
    % number of beam elements for this beam. 
    n_elem = n_elems(i);
    % lambda multipliers for all beam elements. 
    lambdas = repmat(lambda, 1, n_elem+1); 
    % kinematic dof for each beam element. 
    beam_dofs = repmat(dof, 1, n_elem +1);
    % constraints for each beam element. 
    beam_constraints = repmat(n_constraints, 1, n_elem + 1);
    % the fixed dof for each beam element. 
    beam_fixed_dofs = zeros(1, n_elem + 1);
    beam_fixed_dofs(1) = fixed_dof;
    [beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elem);
    beam_directors = repmat(vertcat(d1, d2, d3), 1, n_elem + 1);
    % beam elements and beam object. 
    beam_elems = create_beam_elements(n_elem, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
    beam = Beam(n_elem, beam_elems, dof, n_constraints); 
    % add to the array of the beams. 
    all_beams = [all_beams, beam]; 
end

% Material properties. 
A = 0.2; % in^2 
I = (1/6000); %  in^4 
E = 30e6; % lb/in^2 

% Material property matrix. 
GJ = E*I;
EI2 = E*I; 
EI3 = E*I;
EA = 2.3077e6; 
GA2 = 2.3077e6; 
GA3 = 2.3077e6; 
C = diag([GA2, GA3, EA, EI2, EI3, GJ]); 

% plot axis limits. 
x_lim1 = -L0 - 50;
x_lim2 = L0 + 50; 
y_lim1 = -50; 
y_lim2 = 50; 
z_lim1 = -50; 
z_lim2 = 50;

% vertical forces. 
F1_z = -1.35; 
F2_z = -0.85; 

% position of the forces. 
p1 = [L0, 0, 0]; 
p2 = [52.03, 0, 0];

% arrays for storing results. 
avg_iterations = []; 
Ub_vals = []; 
Vb_vals = []; 
Ub_norms = [];  
Vb_norms = []; 

% analytical Ub, Vb. 
Ub_analytical = 30.75; 
Vb_analytical = 66.96; 
%% SIMULATIONS. 

% iterate over all the beams. 
for i=1:length(all_beams)
    % extract the beam. 
    beam = all_beams(i); 
    dL = L0/beam.n_elements; 

    % node closest to p2. 
    node_p2 = p2(1)/dL; 
    if (i == 1)
        node_p2 = 2; 
    else

        if (node_p2 - floor(node_p2) < 0.5)
            node_p2 = floor(node_p2);
        else
            node_p2 = ceil(node_p2);
        end
    end
    
    % external force vector. 
    f_ext = zeros((beam.n_nodes) * beam.dof_per_node, 1); 
    f_ext(node_p2 * beam.dof_per_node + 3) = F2_z; % force at the middle of the beam. 
    f_ext((beam.n_elements) * beam.dof_per_node + 3) = F1_z; % force at the free end of the beam. 
    
    % Gauss quadrature points. 
    n_gauss_points = 1; 
    
    % max iter count. 
    max_iter = 50; 
    
    % Tolerance. 
    Tol = 1e-7; 
    
    % fixed and free nodes. 
    fixed_nodes = [1]; 
    free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, beam.n_elements + 1);
    
    % load steps. 
    load_steps = 100; 
    total_iterations = 0;  
    for j=1:load_steps
        force = f_ext * (j/load_steps); 
        [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, force, 1, free_dof);
        total_iterations = total_iterations + iter;         
    end
    
    end_beam = beam.beam_elements(beam.n_elements); 
    % intial end beam position. 
    X_initial = end_beam.x2; 
    % updated position. 
    X_current = end_beam.x2_t;  
    % calculate the difference. 
    dX = X_current - X_initial;
    % display end node translation. 
    % disp("Mean iterations: " + num2str(round(avg_iterations, 2)) + " U_b = " + num2str(abs(dX(1))) + ", V_b = " + num2str(abs(dX(3)))); 
    U_b = abs(dX(1)); 
    V_b = abs(dX(3)); 
    Ub_vals = [Ub_vals, U_b]; 
    Vb_vals = [Vb_vals, V_b];
    delta_Ub = norm(U_b - Ub_analytical);  
    delta_Vb = norm(V_b - Vb_analytical); 
    Ub_norms = [Ub_norms, delta_Ub/norm(Ub_analytical)]; 
    Vb_norms = [Vb_norms, delta_Vb/norm(Vb_analytical)]; 
    avg_iter = total_iterations / load_steps; 
    avg_iterations = [avg_iterations, avg_iter];  

end

%% PLOTS. 

% Plot UB vs number of beam elements. 
figure(1); 
plot(n_elems, Ub_vals,  "-o",color="red", LineWidth=1); 
yline(Ub_analytical, color="blue", LineWidth=1); 
xlabel("Number of Beam elements [-]"); 
ylabel("Ub [m]");
title("Horizontal displacement Ub vs No. beam elements"); 
legend(["Ub", "Ub analytical"]); 
grid; 

figure(2); 
loglog(n_elems, Ub_norms,  "-o",color="red", LineWidth=1); 
xlabel("Number of Beam elements [-]"); 
ylabel("||ΔUb|| / ||Ub analytical||");
title("||ΔUb|| / ||Ub analytical|| vs No. Beam elements"); 
legend(["||ΔUb|| / ||Ub analytical||"]); 
grid; 

figure(3); 
plot(n_elems, Vb_vals,  "-o",color="red", LineWidth=1); 
yline(Vb_analytical, color="blue", LineWidth=1); 
xlabel("Number of Beam elements [-]"); 
ylabel("Vb [m]");
title("Vertical displacement Vb vs No. Beam elements"); 
legend(["Vb", "Vb analytical"]); 
grid; 

figure(4); 
loglog(n_elems, Vb_norms,  "-o",color="red", LineWidth=1); 
xlabel("Number of Beam elements [-]"); 
ylabel("||ΔVb|| / ||Vb analytical||");
title("||ΔVb|| / ||Vb analytical|| vs No. Beam elements"); 
legend(["||ΔVb|| / ||Vb analytical||"]); 
grid; 
