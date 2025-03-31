
% returns a vector of dd_i - variations of the directors for a given beam
% node.
function [delta_d] = compute_delta_d(d1, d2, d3, d_theta)
    delta_d1 = -1 * cross(d1, d_theta); 
    delta_d2 = -1 * cross(d2, d_theta); 
    delta_d3 = -1 * cross(d3, d_theta);
    delta_d = [delta_d1; delta_d2; delta_d3]; 
end