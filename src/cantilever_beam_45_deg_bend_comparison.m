% 45 degree bend - cantilever beam, load cases.
clear;
clc;
%% CONSTANTS. 
n_elems = 10; % number of beam elements. 

dof = 12; % kinematic degrees of freedom per node. 
fixed_dof = 12; % fix the kinematic dof at the first node. 
n_constraints = 6; % number of constraints per node. 
lambda = zeros(n_constraints,1); % lambda multipliers per node - initialized as zeros. 

% lambda multipliers for all beam elements. 
lambdas = repmat(lambda, 1, n_elems+1); 
% kinematic dof for each beam element. 
beam_dofs = repmat(dof, 1, n_elems +1);
% constraints for each beam element. 
beam_constraints = repmat(n_constraints, 1, n_elems + 1);
% the fixed dof for each beam element. 
beam_fixed_dofs = zeros(1, n_elems + 1);
beam_fixed_dofs(1) = fixed_dof;

% discretize the domain into beam coordinates, beam lengths and directors. 
% [beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);

% beam elements and beam object. 
all_beams = [];

n_cases = 5; 
for i=1:n_cases
    Radius = 100;
    bend_angle = (45/180) * pi; % 45 degree bend. 
    [beam_coordinates, beam_directors, beam_lengths] = discretize_circle(Radius, bend_angle, n_elems);
    beam_elems = create_curved_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
    beam = Beam(n_elems, beam_elems, dof, n_constraints); 
    all_beams = [all_beams, beam];
end

% Material properties. 
A = 1; 
I = (1/12); 
J = (1/6);
E = 1e7; 

% Material property matrix. 
GJ = E*I;
EI2 = E*I; 
EI3 = E*I;
EA = E*A; 
GA2 = 5e6; 
GA3 = 5e6; 
C = diag([GA2, GA3, EA, EI2, EI3, GJ]);  

% Gauss quadrature points. 
n_gauss_points = 1; 

% max iter count. 
max_iter = 50; 

% fixed and free nodes. 
fixed_nodes = [1]; 
free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, n_elems + 1); 

% displacement coordinates for the beams. 
Beam_displacements = []; 

% Tolerance. 
Tol = 1e-7; 
%% SIMULATION. 
load_steps = 10; 
% force components - end node. 
fx_vals = [0, 0, 0, 0, 0]; 
fy_vals = [0, 0, 0, 0, 0];  
fz_vals = [0, 150, 300, 450, 600]; 

for i=1:n_cases
    % beam object. 
    beam = all_beams(i); 

    % force node positions. 
    node_fx = beam.n_elements * beam.dof_per_node + 1; 
    node_fy = beam.n_elements * beam.dof_per_node + 2;
    node_fz = beam.n_elements * beam.dof_per_node + 3; 
    
    fx = fx_vals(i); 
    fy = fy_vals(i); 
    fz = fz_vals(i); 
    
    output = "Forces: fx = " + num2str(fx) + "N, Fy = " + num2str(fy) + "N, Fz = " + num2str(fz);   
    disp(output); 

    % iterate over each load step. 
    for j=1:load_steps
        f_ext = zeros((beam.n_nodes) * beam.dof_per_node, 1); 
        f_ext(node_fx) = fx * (j/load_steps); 
        f_ext(node_fy) = fy * (j/load_steps); 
        f_ext(node_fz) = fz * (j/load_steps);
        % solve using Newton-Rhapson method. 
        [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, f_ext, 1, free_dof); 
    end
    disp(iter); 
    Beam_coordinates = beam.get_deformed_beam_coords();
    Beam_displacements = [Beam_displacements; Beam_coordinates]; 

end

%% PLOTS. 
colors = ["r", "g", "b", "c", "m"]; 
legend_names = ["0N", "150N", "300N", "450N", "600N"]; 
for i=1:n_cases
    i_start = 1 + (i-1) * 3; 
    i_end = i_start + 2; 
    coords = Beam_displacements(i_start:i_end,:);
    X = coords(1,:); 
    Y = coords(2,:); 
    Z = coords(3,:); 
    plot3(X, Y, Z, color=colors(i), Marker="o", MarkerFaceColor="blue", MarkerSize=2, LineWidth=1); 
    hold on; 
end
grid on; 
legend(legend_names); 
title("Cantilever Beam 45 deg bend - Number of elements = " + num2str(n_elems));
xlabel("x-axis [m]"); 
ylabel("y-axis [m]"); 
zlabel("z-axis [m]"); 
xlim([-10, 60]); 
ylim([-50, 80]); 
zlim([0, 80]); 
