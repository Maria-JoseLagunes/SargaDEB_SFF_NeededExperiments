%% Temperature correction dynamics
% Maria José Lagunes
% First created : 2024/01/08
% Last modified on: 2024/01/08

% Tests different values of parameters to observe what happens if whe
% change one parameter. Dynamics of Arrhenius equation in c(T) which
% evaluates the effect of temperature in metabolic rates

% The 1 parameter version calculates the simple Arrhenius correction
% factor based on different Arrhenius temperatures (T_A) and different
% physiological/environnment temperatures

% The 3 parameters version calculates the correction factor according to
% the upper or lower temperature tolerance paramters. Having 3 parameters
% that could vary : T_A, T_AH, T_H or T_A, T_AL, T_L, respectively.
% Boundary temperature is lower or higher depending on reference
% temperature (T_ref)

% The 5 parameters version calculates the correction factor based on
% temperatures within both temperature boundaries ans has paramters T_A,
% T_AH, T_H, T_AL, T_L

% Input
% T : K, vector of different temperatures
% params_T : K, structure with 1,3 or 5 for the different parameters + reference
% temperature (T_ref)
% arrheniusEquation : scalar which indicates the type of Arrhenius equation being
% evaluated

% Output
% TC correction factor as a function of temperatures that will affect
% physiological rates

function c_T = Arrhenius_equation(temperatureParameters, TC, arrheniusModel)

%% Extract parameters values from structure 

T_A = temperatureParameters.T_A ; % K, Arrhenius temperature
T_ref = temperatureParameters.T_ref; % K, reference temperature



%% Convert temperature tested in K
T = TC + 273.15 ; % K, temperature tested (physiological/environmental)

%% Arrhenius equations according to paramters
if arrheniusModel == "simple"       %Simple Arrhenius equation
    s_A = exp (T_A / T_ref  - T_A ./ T) ; % simple Arrhenius factor
    c_T = s_A ;

elseif arrheniusModel == "upper"     
    T_H = temperatureParameters.T_H; % K, upper boundary tolerance temperature
    T_AH = temperatureParameters.T_AH; % K, upper boundary Arrhenius temperature
    s_A = exp (T_A / T_ref  - T_A ./ T) ; % simple Arrhenius factor
    s_H_up = 1 + exp (  T_AH / T_H   - T_AH / T_ref ); % Upper boundary Arrhenius equation (High temperature), upper part of ratio
    s_H_down= 1 + exp(  T_AH / T_H   - T_AH ./ T) ;% Upper boundary Arrhenius equation(High temperature), down  part of ratio
    s_H_ratio = s_H_up ./ s_H_down  ;
    c_T = s_A .* s_H_ratio ;

elseif arrheniusModel == "lower"     
    T_L = temperatureParameters.T_L; % K, lower boundary tolerance temperature
    T_AL = temperatureParameters.T_AL; % K, lower boundary Arrhenius temperature
    s_A = exp (T_A / T_ref  - T_A ./ T) ; % simple Arrhenius factor
    s_L_up = 1 + exp (  T_AL / T_ref  - T_AL / T_L); % Lower boundary Arrhenius equation (Low temperature), upper part of ratio
    s_L_down = 1 + exp(  T_AL ./ T  - T_AL / T_L); % Lower boundary Arrhenius equation (Low temperature), down part of ratio
    s_L_ratio = s_L_up ./ s_L_down ;    
    c_T = s_A .* s_L_ratio ;

elseif arrheniusModel == "both"         % Extended Arrhenius equation
    T_H = temperatureParameters.T_H; % K, upper boundary tolerance temperature
    T_AH = temperatureParameters.T_AH; % K, upper boundary Arrhenius temperature
    T_L = temperatureParameters.T_L ; % K, lower boundary tolerance temperature
    T_AL = temperatureParameters.T_AL; % K, lower boundary Arrhenius temperature
    % s_A = exp (T_A / T_ref  - T_A ./ T) ; % simple Arrhenius factor
    % s_H_up = 1 + exp (  T_AH / T_H   - T_AH / T_ref ); % Upper boundary Arrhenius equation (High temperature), upper part of ratio
    % s_H_down= 1 + exp(  T_AH / T_H   - T_AH ./ T) ;% Upper boundary Arrhenius equation(High temperature), down  part of ratio
    % s_H_ratio = s_H_up ./ s_H_down  ;
    % s_L_up = 1 + exp (  T_AL / T_ref  - T_AL / T_L); % Lower boundary Arrhenius equation (Low temperature), upper part of ratio
    % s_L_down = 1 + exp(  T_AL ./ T  - T_AL / T_L); % Lower boundary Arrhenius equation (Low temperature), down part of ratio
    % s_L_ratio = s_L_up ./ s_L_down ; 
    % c_T = s_A .* s_L_ratio .* s_H_ratio ; % -, Temperature correction factor
    c_T = exp(T_A/T_ref - T_A./T) * ...
    (1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H  -  T_AH/T_ref)) ./ ...
(1 + exp(T_AL./T - T_AL/T_L) + exp(T_AH/T_H  -  T_AH./T)) ; % -, temperature correction factor


end

end





