
% 45 degree bend - cantilever beam, convergence plot. 
clear;
clc;
%% CONSTANTS. 
n_elems = [5, 10, 20, 40]; % number of beam elements. 

% discretize the domain into beam coordinates, beam lengths and directors. 
% [beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);

% beam elements and beam object. 
all_beams = [];

n_cases = length(n_elems); 
for i=1:n_cases
    Radius = 100;
    dof = 12; % kinematic degrees of freedom per node. 
    fixed_dof = 12; % fix the kinematic dof at the first node. 
    n_constraints = 6; % number of constraints per node. 
    lambda = zeros(n_constraints,1); % lambda multipliers per node - initialized as zeros. 
    bend_angle = (45/180) * pi; % 45 degree bend. 
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
    [beam_coordinates, beam_directors, beam_lengths] = discretize_circle(Radius, bend_angle, n_elem);
    beam_elems = create_curved_beam_elements(n_elem, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
    beam = Beam(n_elem, beam_elems, dof, n_constraints); 
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

% Tolerance. 
Tol = 1e-7; 
%% SIMULATION. 
load_steps = 10; 
indx = 2; 

% force components - end node. 
fx_vals = [0, 0, 0, 0, 0]; 
fy_vals = [0, 0, 0, 0, 0];  
fz_vals = [300, 450, 600]; % load levels. 

tip_coordinates = []; 

% tip coordinates - Cardona. 
cardona_X = [22.14, 18.38, 15.55]; 
cardona_Y = [58.64, 52.11, 47.04]; 
cardona_Z = [40.35, 48.59, 53.50]; 

% tip coordinates ref. [BB79] 
BB79_X = [22.33, 18.62, 15.79]; 
BB79_Y = [58.84, 52.32, 47.23]; 
BB79_Z = [40.08, 48.39, 53.37]; 

% tip coordinate ref. [SVQ86] 
SVQ86_X = [22.5, 0, 15.9]; 
SVQ86_Y = [59.2, 0, 47.2]; 
SVQ86_Z = [39.5, 0, 53.4]; 

for i=1:n_cases
    % beam object. 
    beam = all_beams(i); 

    % fixed and free nodes. 
    fixed_nodes = [1]; 
    free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, beam.n_elements + 1);

    % force node positions. 
    node_fx = beam.n_elements * beam.dof_per_node + 1; 
    node_fy = beam.n_elements * beam.dof_per_node + 2;
    node_fz = beam.n_elements * beam.dof_per_node + 3; 
    
    fx = fx_vals(indx); 
    fy = fy_vals(indx); 
    fz = fz_vals(indx);  
    % fz = fz_vals(i);  

    % iterate over each load step. 
    for j=1:load_steps
        f_ext = zeros((beam.n_nodes) * beam.dof_per_node, 1); 
        f_ext(node_fx) = fx * (j/load_steps); 
        f_ext(node_fy) = fy * (j/load_steps); 
        f_ext(node_fz) = fz * (j/load_steps);
        % solve using Newton-Rhapson method. 
        [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, f_ext, 1, free_dof); 
    end
    
    Beam_coordinates = beam.get_deformed_beam_coords();
    coords = Beam_coordinates(:,beam.n_elements + 1);
    tip_coordinates = [tip_coordinates, coords]; 
end

%% PLOTS. 
colors = ["r", "g", "b", "c", "m"]; 
legend_names = ["Numerical", "Cardona", "SVQ86", "BB79"]; 
legend_names_2 = ["Numerical", "Cardona", "BB79"]; 

% plot X pos as function of number of elements. 
figure(1); 
X_numerical = tip_coordinates(1,:);  
plot(n_elems, X_numerical, color="red", Marker="o",MarkerFaceColor="black", MarkerSize=2, LineWidth=1)
hold on; 
yline(cardona_X(indx),'g', LineWidth=1);
hold on; 
if (indx ~= 2)
    yline(SVQ86_X(indx),'b', LineWidth=1);
end
hold on; 
yline(BB79_X(indx),'c', LineWidth=1);
grid on; 
if (indx ~= 2)
    legend(legend_names); 
else
    legend([legend_names_2]);
end
title("X tip coordinate, Load case = " + num2str(fz_vals(indx)) + "N");
xlabel("beam elements [-]");
ylabel("X [m]");

figure(2);  
Y_numerical = tip_coordinates(2,:);  
plot(n_elems, Y_numerical, color="red", Marker="o",MarkerFaceColor="black", MarkerSize=2, LineWidth=1)
hold on; 
yline(cardona_Y(indx),'g', LineWidth=1);
hold on; 
if (indx ~= 2)
    yline(SVQ86_Y(indx),'b', LineWidth=1);
end
hold on; 
yline(BB79_Y(indx),'c', LineWidth=1);
grid on; 
if (indx ~= 2)
    legend(legend_names, Location="southeast"); 
else
    legend([legend_names_2], Location="southeast");
end
title("Y tip coordinate, Load case = " + num2str(fz_vals(indx)) + "N");
xlabel("beam elements [-]");
ylabel("Y [m]");


figure(3);  
Z_numerical = tip_coordinates(3,:);  
plot(n_elems, Z_numerical, color="red", Marker="o",MarkerFaceColor="black", MarkerSize=2, LineWidth=1)
hold on; 
yline(cardona_Z(indx),'g', LineWidth=1);
hold on; 
if (indx ~= 2)
    yline(SVQ86_Z(indx),'b', LineWidth=1);
end
hold on; 
yline(BB79_Z(indx),'c', LineWidth=1);
grid on; 
if (indx ~= 2)
    legend(legend_names); 
else
    legend([legend_names_2]);
end 
title("Z tip coordinate, Load case = " + num2str(fz_vals(indx)) + "N");
xlabel("beam elements [-]");
ylabel("Z [m]");



