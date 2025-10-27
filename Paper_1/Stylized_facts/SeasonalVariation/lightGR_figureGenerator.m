% Light intenisty, Relative growth rate (doublings day-1)
load("simu.mat")
load("obs.mat")



MGSF3_SPR= [0.0781902552204176
0.080139211136891
0.104083526682135
0.0584222737819026
]; 



MGSF3_SUM = [0.0794431554524362
0.0738747099767981
0.0826450116009281
0.060092807424594
]; 


MG = [MGSF3_SPR , MGSF3_SUM];





figure(5)
for i = 1:length(obs)
    % Convert forcing variables in known units
    I = simu(i).lightIntensity / 1e-6 / 3600 ; % micro mol E m-2 s-1
    lightTested = categorical(I); 
    bar(lightTested, obs(i).r_D, 'FaceColor',simu(i).col); hold on; 

    scatter( lightTested, MG(2,i)) % row 2 25




    
    xlabel("Temperature (°C)");
    ylabel('Growth rate (doublings day^{-1})');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
end



diffLight = simu(2).lightIntensity * 100 / simu(1).lightIntensity; 
diffGrowth =  (obs(2).r_D * 100 / obs(1).r_D) ; 


diffGrowth_MG = MG(2,2) * 100 / MG(2,1) ; 







	