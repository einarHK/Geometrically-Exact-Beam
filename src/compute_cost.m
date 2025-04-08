
% cost function - computes the cost between computed stress-strain pair and
% the stress-strain pair coming from data. 
% INPUT: e - strain computed, s - stress computed, e_data - strain from data, s_data - stress from data, C - weight matrix.   
% OUTPUT: Cost J. 
function J = compute_cost(e, s, e_data, s_data, C)
    delta_e = C*e - C*e_data; 
    delta_s = inv(C) * s - inv(C) * s_data; 
    J = 0.5 * dot(delta_e, delta_e) + 0.5 * dot(delta_s, delta_s); 
end


