
clear; 
clc; 
%% CONSTANTS. 
n_elems = 10; % number of beam elements. 
L0 = 102.75; % initial length of the beam. 
dL = L0 / n_elems; %

x1 = [0, 0, 0]'; % start position of the beam.
x2 = [L0, 0, 0]'; % end position of the beam. 
d3 = [1, 0, 0]'; % axial director. 
d2 = [0, 1, 0]'; % shear directors 2, 1. 
d1 = [0, 0, 1]'; 
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
[beam_coordinates, beam_lengths] = discretize_domain(x1, x2, n_elems);
beam_directors = repmat(vertcat(d1, d2, d3), 1, n_elems + 1);

% beam elements and beam object. 
beam_elems = create_beam_elements(n_elems, beam_coordinates, beam_directors, beam_lengths, beam_constraints, beam_dofs, beam_fixed_dofs, lambdas);
% beam = Beam(n_elems, beam_elems, dof, n_constraints); 

% Material properties. 
A = 0.2; % in^2 
I = (1/6000); % in^4 
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
x_lim = L0 + 50;
y_lim = 50; 
z_lim = 50; 

% vertical forces. 
F1_z = -1.35; 
F2_z = -0.85; 

% position of the forces. 
p1 = [L0, 0, 0]; 
p2 = [52.03, 0, 0];

% node closest to p2. 
node_p2 = p2(1)/dL; 
if (node_p2 - floor(node_p2) < 0.5)
    node_p2 = floor(node_p2);
else
    node_p2 = ceil(node_p2);
end

% Gauss quadrature points. 
n_gauss_points = 1; 

% max iter count. 
max_iter = 50; 

% Tolerance. 
Tol = 1e-7; 

% forces for each case. 
F1_forces = [0.2 * F1_z, 0.4 * F1_z, 0.6 * F1_z, 0.8 * F1_z, F1_z]; % force applied at free end. 
F2_forces = [0.2 * F2_z, 0.4 * F2_z, 0.6 * F2_z, 0.8 * F2_z, F2_z]; % force applied at the middle of the beam. 

% create beam class objects for each force case. 
all_beams = [];  

for i=1:length(F1_forces)
    beam = Beam(n_elems, beam_elems, dof, n_constraints); 
    all_beams = [all_beams, beam]; 
end

% end node - free. 
fixed_end_node = 0; 
%% Run Simulation. 
% beam coordinates. 
X = []; 
Y = []; 
Z = [];
load_steps = 10; % load steps for each force. 
for i=1:length(F1_forces)
    f1_z = F1_forces(i); 
    f2_z = F2_forces(i); 
    beam = all_beams(i); 
    % external force vector. 
    for j=1:load_steps
        f_ext = zeros((beam.n_nodes) * beam.dof_per_node, 1); 
        f_ext(node_p2 * beam.dof_per_node + 3) = f2_z * (j/load_steps); % force at the middle of the beam. 
        f_ext((beam.n_elements) * beam.dof_per_node + 3) = f1_z * (j/load_steps); % force at the free end of the beam. 
        [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, f_ext, 1, fixed_end_node); % iterative solver. 
    end
    [X1, Y1, Z1] = beam.get_beam_coords(); 
    X = [X;X1]; 
    Y = [Y;Y1]; 
    Z = [Z;Z1];
end
%% Plot results. 
colors = ["r", "g", "b", "c", "m"];
legend_names = ["20% P", "40% P", "60% P", "80% P", "100% P"];

for i=1:length(colors)
    x = X(i,:);
    y = Y(i,:); 
    z = Z(i,:); 
    color = colors(i);
    plot(x, z, LineStyle='-', Color=color, Marker="o", MarkerFaceColor="black", LineWidth=1);
    hold on; 
end
legend(legend_names); 
grid on; 
xlabel("x-axis"); 
ylabel("z-axis");
xlim([0, 120]); 
ylim([-80,20]);
title("Successive configurations - Cantilever beam");

