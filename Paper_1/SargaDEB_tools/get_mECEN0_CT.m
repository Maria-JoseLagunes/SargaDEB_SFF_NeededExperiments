function [m_EC_0, m_EN0] = get_mECEN0_CT(rmax, pars_m0)


m_EC_0 = (pars_m0.j_ECAm - pars_m0.kap_EC * pars_m0.j_ECM - pars_m0.kap_EC * pars_m0.y_ECV * rmax)/ ...
                ((1 - pars_m0.kap_EC) * pars_m0.k_EC + pars_m0.kap_EC  * rmax);

m_EN_0 = (pars_m0.j_ENAm - pars_m0.kap_EN * pars_m0.j_ENM - pars_m0.kap_EN * pars_m0.y_ENV * rmax)/ ...
                ((1 - pars_m0.kap_EN) * pars_m0.k_EN + pars_m0.kap_EN * rmax);


end


 
