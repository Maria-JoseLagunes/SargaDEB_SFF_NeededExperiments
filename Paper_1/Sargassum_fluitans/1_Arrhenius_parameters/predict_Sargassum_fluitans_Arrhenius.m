function [prdData, info] = predict_Sargassum_fluitans_Arrhenius(pars, data, auxData)


% data and auxData

vars_pull(data);
vars_pull(auxData);

if pars.T_L > pars.T_ref || pars.T_H < pars.T_ref
  info = 0; prdData = []; return
else
  info = 1;
end



%compute temperature correction factors
pars_T = [pars.T_A; pars.T_L; pars.T_H; pars.T_AL; pars.T_AH]; % K 


TC = tempcorr(C2K(data.TC_RG(:,1)), pars.T_ref, pars_T);  %-

prdData.TC_RG = TC; 



end


