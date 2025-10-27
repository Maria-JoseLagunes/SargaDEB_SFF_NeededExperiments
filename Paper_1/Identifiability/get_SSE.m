% Calculate the Sum of the Squared Error (SSE)
% Input
% theta value : initial parameter/ baseline parameter
% var_value : value to evaluate after estimation from new seed point
% Output
% SSE = Sum of the Squared Error for each free parameter tested



function SSE = get_SSE(theta_value,var_value)
    sign_value = (var_value - theta_value); 
    sign_value = sign(sign_value); 
    SSE = sign_value * (theta_value - var_value)^2 / ((theta_value)^2 + (var_value)^2);
end