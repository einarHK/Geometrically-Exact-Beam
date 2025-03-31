
% displays the stress strain pair for a given beam element. 
function display_stress_strain(stress, strain)
    
    fprintf("Comp:    strain:     stress:\n");
    for i=1:length(stress)
       fprintf("%d      %f     %f\n", i, strain(i), stress(i));
    end
end