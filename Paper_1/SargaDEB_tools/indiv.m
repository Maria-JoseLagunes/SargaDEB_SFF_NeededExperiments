% ---------------------------------------------------------------------
%% Compute numerical integration
% Integration using solver
% Output
% [t, mECENV] = vector of t and EC EN MV
% Maria José Lagunes
% First created : 2024/03/12
% ---------------------------------------------------------------------
function [t, mECENV, J_struct] = indiv(simulationConditions)


forcingVariables.CO2 = simulationConditions.CO2; 
forcingVariables.lightIntensity = simulationConditions.lightIntensity; 
forcingVariables.N = simulationConditions.N; 
forcingVariables.P = simulationConditions.P; 
forcingVariables.photoPeriod = simulationConditions.photoPeriod; 


t0 = simulationConditions.t0; 
tFinal = simulationConditions.tFinal;
ct = simulationConditions.ct; 
initialConditions = simulationConditions.initialConditions; 
timeSimulation = [t0:tFinal]' ;

if  sum(contains( fieldnames( simulationConditions) , 'tSpinTune' )) ==  1
    tSpinTune = simulationConditions.tSpinTune; 
    timeSimulation_SpineTune = [t0:tSpinTune ]' ;
    options = odeset('RelTol',1e-10,'AbsTol',1e-12);
    
    % Integration pour spin-up
    [t, mECENV] = ode89(@(t,mECENV) get_mECENV(t,  mECENV, simulationConditions.par, ct, forcingVariables), timeSimulation_SpineTune, initialConditions, options);
    %% With initial conditons at equilibrium
    initialConditions_equilibre = mECENV(end,:);
    initialConditions_equilibre(1, [1 2]) = mECENV(end,[1 2]); %Retrieving state variables final state, to use as intial condition in simulations
    
    M_V_0 = simulationConditions.dW0 / (simulationConditions.par.w_V +  ...
        (initialConditions_equilibre(1, 2) * simulationConditions.par.w_EN) +  ...
        (initialConditions_equilibre(1, 1) * simulationConditions.par.w_EC)); 
    
    initialConditions_equilibre(1, 3) = M_V_0; %Retrieving state variables final state, to use as intial condition in simulations
    

    [t, mECENV] = ode89(@(t,mECENV) ...
        get_mECENV(t,  mECENV, simulationConditions.par, ct, forcingVariables), ...
        timeSimulation, initialConditions_equilibre, options);


else

    options = odeset('RelTol',1e-10,'AbsTol',1e-12);

    % Integration pour spin-up
    [t, mECENV] = ode89(@(t,mECENV) ...
        get_mECENV(t,  mECENV, simulationConditions.par, ct, forcingVariables), ...
        timeSimulation, initialConditions, options);


end

%Saving flux as datatable
colName = string({"ct", "j_CO2", "j_I", "j_EC_A", "j_EN_A", "j_EC_C", "j_EN_C", "j_EC_MC", "j_EN_MN", "j_V_MC", ...
                    "j_V_MN", "j_V_M", "j_EC_G", "j_EN_G", "j_VG", "j_EC_R", "j_EN_R", "r", "I", "CO2", "N","P"});

J = zeros(length(t),length(colName)); %initializing flux of organic vector

for i = 1:length(t)
    [~,J(i,:)] = get_mECENV (t(i), mECENV(i,:),simulationConditions.par, ct, forcingVariables); 
end


% Create a structure
J_struct= struct();
% 
% Assign each column of J to a field in the structure
for i = 1:length(colName)
    J_struct.(colName(i)) = J(:, i);
end


end


