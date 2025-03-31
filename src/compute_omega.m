
% calculates omega vector for given beam element. 
% INPUT: interpolation value s, initial beam length L0 and directors for
% node 1 denoted as A and node 2, denoted as B. 
% OUTPUT: omega vector.
function omega = compute_omega(s, L0, d1_At, d2_At, d3_At, d1_Bt, d2_Bt, d3_Bt)
    n1 = 0.5 * (1 - s); % shape function N1. 
    n2 = 0.5 * (1 + s); % shape function N2.
    dn1 = -1/L0; % derivative of shape function N1. 
    dn2 = 1/L0; % derivative of shape function N2. 
    d1 = n1 * d1_At + n2 * d1_Bt; % interpolated d1 value. 
    d2 = n1 * d2_At + n2 * d2_Bt; % interpolated d2 value. 
    d3 = n1 * d3_At + n2 * d3_Bt; % interpolated d3 value. 
    delta_d1 = dn1 * d1_At + dn2 * d1_Bt; % interpolated delta d1 value. 
    delta_d2 = dn1 * d2_At + dn2 * d2_Bt; % interpolated delta d2 value. 
    delta_d3 = dn1 * d3_At + dn2 * d3_Bt; % interpolated delta d1 value.
    omega1 = 0.5 * (transpose(d3)*delta_d2 - transpose(d2)*delta_d3); % first component of omega vector. 
    omega2 = 0.5 * (transpose(d1)*delta_d3 - transpose(d3)*delta_d1); % second component of omega vector. 
    omega3 = 0.5 * (transpose(d2)*delta_d1 - transpose(d1)*delta_d2); % third component of omega vector. 
    omega = [omega1; omega2; omega3]; % omega vector. 
end