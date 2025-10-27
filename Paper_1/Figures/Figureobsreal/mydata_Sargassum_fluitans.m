function [data, auxData, metaData, txtData, weights] = mydata_Sargassum_fluitans

%% set metaData
metaData.phylum     = 'Ochrophyta'; 
metaData.class      = 'Phaeophycaeae'; 
metaData.order      = 'Fucales'; 
metaData.family     = 'Sargassaceae';
metaData.species    = 'Sargassum_fluitans'; 
metaData.species_en = 'SFF'; 


metaData.T_typical  = C2K(25); % K, 


% uni-variate data
% MG2023a
%time (days), total wet weight, temperature, 
% data_t_Ww_28 = [ ... % d, g Ww, °C
%                  0     12     29    ;
%                  5     14     28.5  ;
%                 10     17     28 ];

data_t_Ww_28 = [ ... % d, g Ww, °C
                 0     12     29    ;
                 2     14     28.5  ;
                 4     15     28 ;
                 6     17     29;
                 8     17.5   28;
                 10    18     28.5
                 ];
      

data.tWw_MG2023a_28 = data_t_Ww_28(1:end, [1 2]); %Select 
units.tWw_MG2023a_28 = {'d', 'g Ww'}; label.tWw_MG2023a_28 = {'days' , 'wet weight'}; 
Ww0.tWw_MG2023a_28 = data_t_Ww_28(1,2); units.Ww0.tWw_MG2023a_28 = 'g Ww'; label.Ww0.tWw_MG2023a_28 = 'Initial wet weight'; 
time.tWw_MG2023a_28 = 10 ;  units.time.tWw_MG2023a_28 = 'd'; label.time.tWw_MG2023a_28 = 'days';


temp.tWw_MG2023a_28 = C2K(data_t_Ww_28(1:end,3)); units.temp.tWw_MG2023a_28 = 'K'; label.temp.tWw_MG2023a_28 = 'temperature';
light.tWw_MG2023a_28 = repelem(mean([435 581]) * 1e-6 * 3600, length(temp.tWw_MG2023a_28))'; units.light.tWw_MG2023a_28 = 'mol gamma m-2 s-1'; label.light.tWw_MG2023a_28 = "light intensity"; 
nitrogen.tWw_MG2023a_28 = repelem((0.06 + 3.2) * 1e-6, length(temp.tWw_MG2023a_28))' ; units.nitrogen.tWw_MG2023a_28 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_MG2023a_28 = "N concentration"; 
CO2.tWw_MG2023a_28 = 0.002 ; units.CO2.tWw_MG2023a_28 = 'mol CO2 L-1'; label.CO2.tWw_MG2023a_28 = "C concentration"; 
phosphorous.tWw_MG2023a_28  =  repelem(0.5 * 1e-6, length(temp.tWw_MG2023a_28))'  ; units.phosphorous.tWw_MG2023a_28 = 'mol P043- L-1 '; label.phosphorous.tWw_MG2023a_28 = "P concentration"; 
photoPeriod.tWw_MG2023a_28   = 1; units.photoPeriod.tWw_MG2023a_28 = '-'; label.photoPeriod.tWw_MG2023a_28 = "Photoperiod"; 



bibkey.tWw_MG2023a_28 = 'magana-gallegos_growth_2023'; 


%%  set weights for all real data
weights = setweights(data, []); %% 

weights.tWw_MG2023a_28 =  weights.tWw_MG2023a_28 * 1;



%% pack auxData and txtData for output
auxData.temp = temp;
auxData.light = light;
auxData.nitrogen = nitrogen;
auxData.CO2 = CO2;
auxData.phosphorous = phosphorous;
auxData.photoPeriod = photoPeriod;


auxData.Ww0.tWw_MG2023a_28 = Ww0.tWw_MG2023a_28;





auxData.time.tWw_MG2023a_28 = time.tWw_MG2023a_28;






txtData.units = units; 
txtData.label = label;
txtData.bibkey = bibkey;


end