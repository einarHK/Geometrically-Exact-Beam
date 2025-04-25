% Quarter of a circular arc with total arc length of 1 m. 
clear; 
clc;
%% CONSTANTS. 
% initial beam length. 
L = 1; 
% number of elements. 
n_elems = 20; 
% kinematic degrees of freedom per node. 
dof = 12; 
% constraints per node. 
n_constraints = 6; 
% Initial lagrange multiplers for node. 
lambda = zeros(n_constraints, 1); 

% lambda multipliers for all beam elements. 
lambdas = repmat(lambda, 1, n_elems+1); 
% kinematic dof for each beam element. 
beam_dofs = repmat(dof, 1, n_elems +1);
% constraints for each beam element. 
beam_constraints = repmat(n_constraints, 1, n_elems + 1);
% the fixed dof for each beam element. 
beam_fixed_dofs = zeros(1, n_elems + 1);
% fix first and last node of the beam. 
fixed_dof = 12;
beam_fixed_dofs(1) = fixed_dof;
beam_fixed_dofs(n_elems + 1) = fixed_dof;  

% material property matrix. 
C = zeros(6); 
C(1,1) = 75; 
C(2,2) = 75; 
C(3,3) = 100; 
C(4,4) = 100; 
C(5,5) = 100; 
C(6,6) = 200; 

% discretize the domain into beam coordinates, beam lengths and directors. 
% [beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
Radius = 2/pi;
bend_angle = (1/2) * pi; % 90 degree bend. 
[beam_coordinates, beam_directors, beam_lengths] = discretize_circle(Radius, bend_angle, n_elems);

% beam elements and beam object. 
beam_elems = create_curved_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
beam = Beam(n_elems, beam_elems, dof, n_constraints); 

% Gauss quadrature points. 
n_gauss_points = 1; 

% max iter count. 
max_iter = 50; 

% Tolerance. 
Tol = 1e-12; 

% dimension for e. 
e_dim = 6; 

% the fixed and free nodes. 
fixed_nodes = [1, n_elems + 1]; 
free_dof = compute_free_dof(fixed_nodes, dof, n_constraints, n_elems + 1); 
%% SIMULATION. 
load_steps = 20; 

% force vectors. 
f1 = transpose([-10, 0, -20]) ; 
f2 = transpose([7.5, -7.5, 15]); 
f3 = transpose([0, 10, -20]); 

% vertical load case. 

% force node positions. 
p2x = beam.dof_per_node + 1; 
p2y = beam.dof_per_node + 2; 
p2z = beam.dof_per_node + 3; 
p3x = beam.dof_per_node * 2 + 1; 
p3y = beam.dof_per_node * 2 + 2; 
p3z = beam.dof_per_node * 2 + 3; 
p4x = beam.dof_per_node * 3 + 1; 
p4y = beam.dof_per_node * 3 + 2; 
p4z = beam.dof_per_node * 3 + 3; 
p8x = beam.dof_per_node * 7 + 1; 
p8y = beam.dof_per_node * 7 + 2; 
p8z = beam.dof_per_node * 7 + 3; 
p9x = beam.dof_per_node * 8 + 1;  
p9y = beam.dof_per_node * 8 + 2; 
p9z = beam.dof_per_node * 8 + 3; 
p10x = beam.dof_per_node * 9 + 1; 
p10y = beam.dof_per_node * 9 + 2; 
p10z = beam.dof_per_node * 9 + 3; 
p11x = beam.dof_per_node * 10 + 1; 
p11y = beam.dof_per_node * 10 + 2; 
p11z = beam.dof_per_node * 10 + 3; 
p12x = beam.dof_per_node * 11 + 1; 
p12y = beam.dof_per_node * 11 + 2; 
p12z = beam.dof_per_node * 11 + 3; 
p13x = beam.dof_per_node * 12 + 1; 
p13y = beam.dof_per_node * 12 + 2; 
p13z = beam.dof_per_node * 12 + 3; 
p14x = beam.dof_per_node * 13 + 1; 
p14y = beam.dof_per_node * 13 + 2; 
p14z = beam.dof_per_node * 13 + 3; 
p18x = beam.dof_per_node * 17 + 1; 
p18y = beam.dof_per_node * 17 + 2; 
p18z = beam.dof_per_node * 17 + 3; 
p19x = beam.dof_per_node * 18 + 1; 
p19y = beam.dof_per_node * 18 + 2; 
p19z = beam.dof_per_node * 18 + 3; 
p20x = beam.dof_per_node * 19 + 1;
p20y = beam.dof_per_node * 19 + 2; 
p20z = beam.dof_per_node * 19 + 3; 

% external force vector. 
f_ext = zeros((beam.n_nodes) * beam.dof_per_node, 1); 

% add forces to correct positions. 
f_ext(p2x:p2z) = f_ext(p2x:p2z) + f1;
f_ext(p3x:p3z) = f_ext(p3x:p3z) + f1; 
f_ext(p4x:p4z) = f_ext(p4x:p4z) + f1; 
f_ext(p8x:p8z) = f_ext(p8x:p8z) + f2; 
f_ext(p9x:p9z) = f_ext(p9x:p9z) + f2; 
f_ext(p10x:p10z) = f_ext(p10x:p10z) + f2; 
f_ext(p11x:p11z) = f_ext(p11x:p11z) + f2; 
f_ext(p12x:p12z) = f_ext(p12x:p12z) + f2; 
f_ext(p13x:p13z) = f_ext(p13x:p13z) + f2; 
f_ext(p14x:p14z) = f_ext(p14x:p14z) + f2; 
f_ext(p18x:p18z) = f_ext(p18x:p18z) + f3; 
f_ext(p19x:p19z) = f_ext(p19x:p19z) + f3; 
f_ext(p20x:p20z) = f_ext(p20x:p20z) + f3; 

beam.display_end_node_pos();
% iterate over each load step. 
for i=1:load_steps
    force = f_ext * (i/load_steps); 
    % solve using Newton-Rhapson method. 
    [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, force, 1, free_dof); 
end
%% PLOT RESULTS. 
% analytical solution - node 11. 
x1_analytical = 0.30620366;
y1_analytical = 0.33041646;
z1_analytical = 0.21515427; 

% analytical d1 director. 
d11_analytical = 0.00761486;
d12_analytical = -0.00761485; 
d13_analytical = 0.99994201; 

% analytical d2 director
d21_analytical = 0.70706900;
d22_analytical = -0.70706899; 
d23_analytical = -0.01076908;

% analytical d3 director. 
d31_analytical = 0.70710999;
d32_analytical = 0.70711000; 
d33_analytical = 0.00000000; 

fprintf(1, '\n');

disp("Analytical solution - Node 11: ");
disp("Postition = (" + num2str(x1_analytical) + "m, " + num2str(y1_analytical) + "m, " + num2str(z1_analytical) + "m )"); 
disp("d1 = (" + num2str(d11_analytical) + ", " + num2str(d12_analytical) + ", " + num2str(d13_analytical) + " )"); 
disp("d2 = (" + num2str(d21_analytical) + ", " + num2str(d22_analytical) + ", " + num2str(d23_analytical) + " )"); 
disp("d3 = (" + num2str(d31_analytical) + ", " + num2str(d32_analytical) + ", " + num2str(d33_analytical) + " )"); 

fprintf(1, '\n');

scale = 0.01; 
x_lim1 = -0.1; 
x_lim2 = 1; 
y_lim1 = -0.1; 
y_lim2 = 0.7; 
z_lim1 = 0; 
z_lim2 = 0.4;

beam.plot_undeformed_state(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, "");
hold on; 
beam.plot_deformed_beam(x_lim1, x_lim2, y_lim1, y_lim2, z_lim1, z_lim2, scale, "Deformed (\color{red}red\color{black}) & Undeformed (\color{blue}blue\color{black})");

node_indices = [11]; 

% compute stress strain pair for the eight element. 
elements = [8]; 
[stress, strain] = beam.compute_stress_strain(elements, C); 
fprintf("Strain and stresses - Beam element %d\n", elements(1));
display_stress_strain(stress, strain); 
fprintf("\n");

deformed_coords = beam.get_deformed_beam_coords(); 
node_coords = deformed_coords(:,node_indices(1)); 
x1 = node_coords(1); 
y1 = node_coords(2); 
z1 = node_coords(3); 
 
dx = abs(x1 - x1_analytical) / abs(x1_analytical); 
dy = abs(y1 - y1_analytical) / abs(y1_analytical); 
dz = abs(z1 - z1_analytical) / abs(z1_analytical);

beam.display_node_pos(node_indices);
beam.display_node_directors(node_indices);

fprintf("\n"); 
fprintf("Relative error - Node position: ");
fprintf("|dx|/|x| = %f , |dy|/|y| = %f , |dz|/|z| = %f \n", dx, dy, dz); 

%% PLOT - Stress component for each beam. 
elements = 1:n_elems; 
[s, e] = beam.compute_stress_strain(elements, C); 
s1 = []; 
s2 = []; 
s3 = []; 
s4 = []; 
s5 = []; 
s6 = []; 

for i=1:n_elems
    i_start = (i - 1) * e_dim + 1;
    i_end = i_start - 1 + e_dim; 
    s_val = s(i_start:i_end); 
    s1 = [s1, s_val(1)]; 
    s2 = [s2, s_val(2)]; 
    s3 = [s3, s_val(3)]; 
    s4 = [s4, s_val(4)]; 
    s5 = [s5, s_val(5)]; 
    s6 = [s6, s_val(6)]; 
end

x_values = 1:n_elems;
stress_components = [s1; s2; s3; s4; s5; s6];
for i=1:6
    figure(i+1); 
    for j=1:n_elems
        s_val = stress_components(i,j);
        x_val = x_values(j); 
        plot([x_val, x_val], [0, s_val], color="blue", LineWidth=1); 
        hold on; 
        scatter(x_val, s_val, MarkerEdgeColor="blue");
        hold on;
    end
    xlabel("Beam element [#]"); 
    if (i < 4)
        ylabel("Stress resultant component " + num2str(i) + " [N/m]");
    else
        ylabel("Stress resultant component " + num2str(i) + " [N]");
    end
    xlim([-1, n_elems + 2]);
    hold off; 
    grid on; 
end



