finalWd_EC = obs.Wd_EC(end); 
finalWd_EN = obs.Wd_EN(end);
finalWd_V = obs.Wd_V(end); 


figure()
final_Wd = [finalWd_EC, finalWd_EN, finalWd_V]; 
names = ["Wd_EC","Wd_EN","Wd_V"];
piechart(final_Wd,names)


finalWw_EC = obs.Wd_EC(end) / dwratio; 
finalWw_EN = obs.Wd_EN(end) / dwratio;
finalWw_V = obs.Wd_V(end) / dwratio; 

initialWw_v = obs.Wd_V(1) / dwratio; 


figure()
final_Ww = [finalWw_EC, finalWw_EN, finalWw_V]; 
names = ["Ww_EC","Ww_EN","Ww_V"];
piechart(final_Ww,names)


% Calculate growth rate using the structure weight
%% Growth rate measure
r_GV = log(finalWw_V /initialWw_v) / obs.timeDays; % d-1, relative growth rate of structure
r_DV = log2(finalWw_V /initialWw_v) / obs.timeDays; % doubling d-1, specific growth rate of the structure
t_DV = 1 / r_DV; % time to double biomass structure

r_GM_V = log(M_V(end) / M_V(1)) / obs.timeDays; % d-1, relative growth rate of structure
% So r_GV and r_GMV is the same. 


r_GM_V_hours = 1/M_V(end)  * (M_V(end)- M_V(end-1)) / (t(end)- t(end-1)); % h-1, relative growth rate of structure


r_GV_hours = log(finalWw_V(end)/(obs.Wd_V(end-1)/dwratio)) / (t(end)- t(end-1)); % h-1, relative growth rate of structure

r = r_GM_V_hours; %in hours

r = (0.1 / 24) * 0.95;  %in hours, rmax from studies


m_EC_0 = (par.j_ECAm - par.kap_EC * par.j_ECM - par.kap_EC * par.y_ECV * r)/ ...
                ((1 - par.kap_EC) * par.k_EC + par.kap_EC  * r);

m_EN_0 = (par.j_ENAm - par.kap_EN * par.j_ENM - par.kap_EN * par.y_ENV * r)/ ...
                ((1 - par.kap_EN) * par.k_EN + par.kap_EN * r);


