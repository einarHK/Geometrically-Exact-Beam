% calculates gamma vector for given beam element. 
% INPUT: interpolation value s, initial beam length L0 and directors for
% node 1 denoted as A and node 2, denoted as B. 
% OUTPUT: omega vector.
function gamma = compute_gamma(s, L0, x1, x2, d1_At, d2_At, d3_At, d1_Bt, d2_Bt, d3_Bt)
    n1 = 0.5 * (1 - s); % n1 shape function. 
    n2 = 0.5 * (1 + s); % n2 shape function. 
    dn1 = -1/L0; % dn1/ds value. 
    dn2 = 1/L0; % dn2/ds value. 
    dx = dn1 * x1 + dn2 * x2; % the derivative of x w.r.t beam axis. 
    d1 = n1 * d1_At + n2 * d1_Bt; % interpolated d1 vector. 
    d2 = n1 * d2_At + n2 * d2_Bt; % interpolated d2 vector. 
    d3 = n1 * d3_At + n2 * d3_Bt; % interpolated d3 vector. 
    gamma1 = transpose(d1) * dx; % gamma 1 value.  
    gamma2 = transpose(d2) * dx; % gamma 2 value. 
    gamma3 = transpose(d3) * dx; % gamma 3 value. 
    gamma = [gamma1; gamma2; gamma3]; % gamma vector.
end