
% function for discretetizing a bent curve. 
function [beam_coords, beam_directors, beam_lengths] = discretize_circle(radius, theta, n_elements)
    beam_directors = [];
    beam_coords = [];
    beam_lengths = [];
    dtheta = theta / n_elements; 
    theta0 = pi; 
    % assume curve lies in the xy plane. 
    for i=1:n_elements + 1
        theta_val = theta0 - (i - 1) * dtheta;
        x = radius * cos(theta_val) + radius; 
        y = radius * sin(theta_val); 
        z = 0; 
        coord = transpose([x, y, z]); 
        beam_coords = [beam_coords, coord]; 
    end
    
    % calculate beam lengths and directors. 
    for i=1:n_elements
        x1 = beam_coords(:,i);
        x2 = beam_coords(:,i+1);
        dx = x2 - x1;
        beam_length = sqrt(transpose(dx) * dx); 
        if (i == 1)
            d1 = [0, 0, 1]; 
            d2 = [1, 0, 0];
            d3 = [0, 1, 0];
            % d2 = [cos(dtheta * (i - 1)), -sin(dtheta * (i - 1)), 0]; % d2 = [1, 0, 0]; % [cos(theta), sin(theta), 0]
            % d3 = [sin(dtheta * (i - 1)), cos(dtheta * (i - 1)), 0]; % d3 = [0, 1, 0]; % [-sin(theta), cos(theta), 0]
            % d3 = cross(d1, d2); 
            d_vector = [transpose(d1); transpose(d2); transpose(d3)]; 
        else
            d1 = [0, 0, 1];
            d3 = dx/beam_length; % d3 director. 
            % d2 = [cos(i * dtheta), -sin(i * dtheta), 0]; 
            % d3 = [sin(dtheta * (i-1)), cos(dtheta * (i-1)), 0]; 
            d2 = cross(d3, d1); % d2 director.
            %d3 = cross(d1, d2); 
            d_vector = [transpose(d1); transpose(d2); d3]; 
            % d_vector = [transpose(d1); transpose(d2); transpose(d3)]; 
        end
        beam_directors = [beam_directors, d_vector]; 
        beam_lengths = [beam_lengths, beam_length];        
    end
    % add director at end node. 
    d1 = [0, 0, 1]; 
    d2 = [cos(theta), -sin(theta), 0]; 
    d3 = [sin(theta), cos(theta), 0]; 
    % d3 = cross(d1, d2);  
    % d2 = cross(d3, d1); 
    d_vector = [transpose(d1); transpose(d2); transpose(d3)];
    beam_directors = [beam_directors, d_vector]; 
end