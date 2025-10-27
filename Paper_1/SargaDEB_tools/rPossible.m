%%
[pars, metaPar, txtPar] = pars_init_Sargassum_fluitans([]);

r_list = linspace(0, 0.5, 500);  
m_EC_0 = zeros(size(r_list));
m_EN_0 = zeros(size(r_list));


for i=1:length(r_list)
    r = r_list(i); 
    m_EC_0(i) = (pars.j_ECAm - pars.kap_EC * pars.j_ECM - pars.kap_EC * pars.y_ECV * r)/ ...
                    ((1 - pars.kap_EC) * pars.k_EC + pars.kap_EC  * r);
    
    m_EN_0(i) = (pars.j_ENAm - pars.kap_EN * pars.j_ENM - pars.kap_EN * pars.y_ENV * r)/ ...
                    ((1 - pars.kap_EN) * pars.k_EN + pars.kap_EN * r);

end

%%
figure(1);
subplot(2,1,1)
plot(r_list, m_EN_0, 'r-', 'LineWidth', 2); hold on;
ylabel('m_{EN0}');
yline(0, 'k--');

subplot(2,1,2)
plot(r_list, m_EC_0, 'b--', 'LineWidth', 2); hold on;
ylabel('m_{EC0}');
yline(0, 'k--');

xlabel('specific growth rate r');
