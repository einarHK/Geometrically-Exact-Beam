
% Cantilever beam, two transvers loads applied at free end and at the middle of
% the cantilever. 
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
beam = Beam(n_elems, beam_elems, dof, n_constraints); 

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
%% Simulation. 
% number of load steps. 
load_steps = 10; 
for i=1:load_steps
    force = (i/load_steps) * f_ext; 
    [iter] = Newtons_method_beam(beam, n_gauss_points, C, max_iter, Tol, force, 1);
end
%% Plot beam. 
beam.plot(x_lim, y_lim, z_lim, 0, "");
% show end node positions and directors. 
beam.display_end_node_directors(); 
beam.display_end_node_pos();
% displacement of the free end node. 
% end beam object.
end_beam = beam.beam_elements(beam.n_elements); 
% intial end beam position. 
X_initial = end_beam.x2; 
% updated position. 
X_current = end_beam.x2_t; 
% calculate the difference. 
dX = X_current - X_initial;
% display end node translation. 
output = "U_b = " + num2str(abs(dX(1))) + ", V_b = " + num2str(abs(dX(3))); 
% compare to the analytical solution. 
Ub_analytical = 30.75; 
Vb_analytical = 66.96; 
output2 = "Ub_analytical = " + num2str(Ub_analytical) + ", Vb_analytical = " + num2str(Vb_analytical); 
% print output. 
display(output); 
display(output2);
