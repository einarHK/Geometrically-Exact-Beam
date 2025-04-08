
% function for computing linearly synthetic (e,s) data. 
% INPUT: Number of data points N_points, C coefficient s=C*e, minimum
% strain e_min, maximum strain e_max, data noise percentage noise. 
% OUTPUT: strain stress vectors e and s. 
function [e,s] = compute_linear_stress_strain(N_points, C, e_min, e_max, noise)
    e = linspace(e_min, e_max, N_points); % strain values. 
    s_vals = C * e; % no noise added. 
    s = zeros(N_points, 1); 
    for i=1:N_points
        s_val = s_vals(i); 
        if (rand(1) < 0.5)
            s_new = s_val - noise * s_val * rand(1); 
            s(i) = s_new; 
        else
            s_new = s_val + noise * s_val * rand(1); 
            s(i) = s_new; 
        end
    end
end
