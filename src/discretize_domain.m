
% function for discretizing a domain into n coordinates and beam lengths. 
function [coordinates, beam_lengths] = discretize_domain(x_start, x_end, n_elems)
    n_points = n_elems + 1; 
    coordinates = [];
    beam_lengths = [];
    dx = (x_end - x_start)/n_elems; 
    for i=1:n_points
        x = x_start + (i-1) * dx; 
        coordinates = [coordinates, x];
        beam_length = sqrt(dx' * dx); 
        beam_lengths = [beam_lengths, beam_length];
    end
end
