
% function for plotting the single beam. 

function plot_beam(beam)
    [X, Y, Z] = beam.compute_coordinates(); % returns node 1 and node 2 X,Y,Z vals. 
    plot3(X, Y, Z, LineStyle='-', Color="red", Marker="o", MarkerFaceColor="blue", LineWidth=1);
    grid on;
    xlabel("x-axis"); 
    ylabel("y-axis");
    zlabel("z-axis");
    xlim([-5, 5]); 
    ylim([-5, 5]); 
    zlim([-5, 5]);
    
end

