%% Description
% Checks if parameter values are in the allowable part of the parameter
%    space of Sargassum DEB model 
% Meant to be run in the estimation procedure
%
% Input
%
%  p: structure with parameters 
%  
% Output
% filter: 0 for hold, 1 for pass
% flag: indicator of reason for not passing the filter
% 1: parameters are not all positive
% 2: one of the kappa and probability parameters is larger than 1
%
% Created by : Maria José Lagunes
% Date modified : 04/12/2024
%%


function [filter, flag] = filter_algal(p)

filter = 0; flag = 0; % default setting of filter and flag

parvec = [p.j_ENAm; p.K_N; p.j_Cm; p.K_C; p.rho_PSU; p.eps_I; p.k_I; p.y_IEC; p.y_CEC; p.y_OI; ...
    p.j_ECAm; p.k_EN; p.k_EC; p.j_ENM; p.j_ECM; p.y_ENV; p.y_ECV; p.kap_EC; p.kap_EN; p.T_A; p.T_AH; p.T_AL];

% parvec = [p.j_ENAm; p.K_N; p.j_CO2m; p.K_C; p.rho_PSU; p.beta_I; p.alpha_I; p.k_I; p.y_IC; p.y_CO2C; p.y_LO2; ...
    % p.j_ECAm; p.k_EN; p.k_EC; p.j_ENM; p.j_ECM; p.y_ENV; p.y_ECV; p.kap_EC; p.kap_EN; p.T_A; p.T_AH; p.T_AL];


 if sum(parvec <= 0) > 0 % all pars must be positive
    flag = 1;
 return;
 end


 probapars = [p.kap_EC; p.kap_EN];
    
    if any(probapars < 0 | probapars > 1)
     flag = 2;
    return;
    end


    if   p.T_H < p.T_ref && p.T_H < p.T_L && p.T_L > p.T_ref
        flag = 3;
    return;
    end


    if p.y_ENV < p.n_NV ||   p.y_ECV < p.n_CV
        flag = 4;
    return;
    end

    if p.j_ECAm < p.j_ECM || p.j_ENAm < p.j_ENM 
        flag = 5;
    return;
    end



    fields = fieldnames(p);
    mEC_fields = fields(startsWith(fields, 'm_EC'));
    mEC_values = cellfun(@(f) p.(f), mEC_fields);
    
    mEN_fields = fields(startsWith(fields, 'm_EN'));
    mEN_values = cellfun(@(f) p.(f), mEN_fields);


    if  any(mEC_values < 0) ||  any(mEN_values < 0)
        flag = 5;
    return;
    end
 
 
filter = 1; 
 
 
end
 
 
 
 
 
 
 % 
% filter = 1;
%  = [p.; p.kap_X; p.kap_P; p.v; p.kap; p.p_M; p.E_G; p.k_J; p.E_Hb; p.E_Hp; p.kap_R; p.h_a; p.T_A];
% parnm = fieldnames(p);
% if sum(contains(parnm,'free')) == 1
%     parnm = fieldnames(p.free);
%         np = numel(parnm);
%     n_par = sum(cell2mat(struct2cell(p.free)));
%     if (n_par == 0)
%     return; % no parameters to iterate
%     end
%     index = 1:np;
%     index = index(cell2mat(struct2cell(p.free)) == 1);  % indices of free parameters
% 
% 
%     q = rmfield(p, 'free'); % copy input parameter matrix into output
%     qvec = cell2mat(struct2cell(q));
%     if sum(qvec <= 0) > 0 % all pars must be positive
%         flag = 1;
%         return; 
%     end
% elseif sum(p <= 0) > 0 % all pars must be positive
%     flag = 1;
%     return;
% end


% %%
% filter = 0; flag = 0; % default setting of filter and flag
% 
% parvec = [p.z; p.kap_X; p.kap_P; p.v; p.kap; p.p_M; p.E_G; p.k_J; p.E_Hb; p.E_Hp; p.kap_R; p.h_a; p.T_A];
% 
%