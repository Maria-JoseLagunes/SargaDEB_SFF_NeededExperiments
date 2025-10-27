% Perform experience based forcing conditons and initial conditions set in
% init file
function [simu, obs] = run_simulation()
    simu = init; 
    
    %Calculate state variables
    for i = 1 : length(simu)
    
        [t, mECENV, J_struct] = indiv(simu(i));
        simu(i).t = t; 
        simu(i).mECENV = mECENV; 
        simu(i).J_struct = J_struct; 
    
    
        % Observables
        obs(i) = get_obs(mECENV,J_struct, simu(i));
        
             
        % Plots
        
       get_plots_2(t,mECENV,J_struct,obs(i), simu(i));
    
    
    end 
end
