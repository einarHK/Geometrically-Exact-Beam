
% Computes mean square error.

function MSE = MSE_error(vector)

MSE = 0; 
for i=1:length(vector)
    MSE = MSE + vector(i) * vector(i);
end
MSE = sqrt(MSE); 

end

