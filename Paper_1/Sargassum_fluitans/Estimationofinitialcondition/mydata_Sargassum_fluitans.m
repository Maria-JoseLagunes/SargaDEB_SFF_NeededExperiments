function [data, auxData, metaData, txtData, weights] = mydata_Sargassum_fluitans

%% set metaData
metaData.phylum     = 'Ochrophyta'; 
metaData.class      = 'Phaeophycaeae'; 
metaData.order      = 'Fucales'; 
metaData.family     = 'Sargassaceae';
metaData.species    = 'Sargassum_fluitans'; 
metaData.species_en = 'SFF'; 


metaData.T_typical  = C2K(25); % K, 

metaData.data_1     = {'T-RGR'; 't-Ww'; 'I-J_O2'}; 

%% set data
% %% Import from tempcorr factor
% load("data_TC_GR.mat"); 
% y
% %% Changing values by %  
% vectorData(1:4,2) =  vectorData(1:4,2) + (vectorData(1:4,2) * 0.2); % Modify value at 28°C to be higher
% data_vector_GR = vectorData; 
% 
% %% Creating data structure
% 
% data.TC_RG = [data_vector_GR(:,1) data_vector_GR(:,2) ];
% units.TC_RG = {'T', '-'}; label.TC_RG = {'°C', 'ct'};
% 
% temp.TC_RG  = C2K(data_vector_GR(:,1)); units.temp.TC_RG = 'K'; label.temp.TC_RG  = 'temperature'; 
% 
% 
% bibkey.TC_RG = 'TC_RG'; 



%% zero-variate data


data.HT2023_CNratio = 32.4;    units.HT2023_CNratio = 'mol C mol N-1';    
label.HT2023_CNratio = 'CN ratio';           bibkey.HT2023_CNratio = 'Hatt2024';
comment.HT2023_CNratio = 'Average from Hatt et al 2024 Figure 4, The green dashed line represents the average nutrient content value for sargasso reported by Lapointe et al. 2021'; 


data.HT2023_CPratio = 886;    units.HT2023_CPratio = 'mol C mol P-1';    
label.HT2023_CPratio = 'CP ratio';           bibkey.HT2023_CPratio  = 'Hatt2024';
comment.HT2023_CPratio = 'Average from Hatt et al 2024 Figure 4, The green dashed line represents the average nutrient content value for sargasso reported by Lapointe et al. 2021'; 


data.HT2023_NPratio = 27.8;    units.HT2023_NPratio = 'mol N mol P-1';    
label.HT2023_NPratio = 'NP ratio';           bibkey.HT2023_NPratio  = 'Hatt2024';
comment.HT2023_NPratio = 'Average from Hatt et al 2024 Figure 4, The green dashed line represents the average nutrient content value for sargasso reported by Lapointe et al. 2021'; 


Ww0.HT2023_CNratio = 2.5; units.Ww0.HT2023_CNratio = 'g Ww'; label.Ww0.HT2023_CNratio = 'Initial wet weight'; 
time.HT2023_CNratio = 500 ;  units.time.HT2023_CNratio = 'd'; label.time.HT2023_CNratio = 'days';


temp.HT2023_CNratio = C2K(25); units.temp.HT2023_CNratio = 'K'; label.temp.HT2023_CNratio = 'temperature';
light.HT2023_CNratio = 500  * 1e-6 * 3600; units.light.HT2023_CNratio = 'mol gamma m-2 s-1'; label.light.HT2023_CNratio = "light intensity"; 
nitrogen.HT2023_CNratio = 1.5 * 1e-6 ; units.nitrogen.HT2023_CNratio = 'mol NO3 + NH4 L-1'; label.nitrogen.HT2023_CNratio = "N concentration"; 
CO2.HT2023_CNratio = 0.002 ; units.CO2.HT2023_CNratio = 'mol CO2 L-1'; label.CO2.tWw_MG2023a_28 = "C concentration"; 
phosphorous.HT2023_CNratio  =  0.10 * 1e-6; units.phosphorous.HT2023_CNratio = 'mol P043- L-1 '; label.phosphorous.HT2023_CNratio = "P concentration"; 
photoPeriod.HT2023_CNratio   = 0; units.photoPeriod.HT2023_CNratio = '-'; label.photoPeriod.HT2023_CNratio = "Photoperiod"; 




% uni-variate data
% HS1987 % DATA FOR SARGASSUM NATANS ! 
data_T_RGR_HS1987 = [ ... %°C ,   doubling d-1
                    12	0.0013589364844904;
                    18	0.0409748892171344;
                    24	0.0480206794682423;
                    30	0.0417725258493353]; 

%12°C
Ww0.tWw_HS1987_12 = 26; units.Ww0.tWw_HS1987_12   = 'g Ww'; label.Ww0.tWw_HS1987_12   = 'Initial wet weight';
time.tWw_HS1987_12  = 21;  units.time.tWw_HS1987_12  = 'd'; label.time.tWw_HS1987_12   = 'days';
WWf_HS1987_12 = 2.^(data_T_RGR_HS1987(1,2) * time.tWw_HS1987_12 ) * Ww0.tWw_HS1987_12  ;
data.tWw_HS1987_12 = [0 Ww0.tWw_HS1987_12 ; time.tWw_HS1987_12   WWf_HS1987_12 ];
units.tWw_HS1987_12 = {'d', 'g Ww'}; label.tWw_HS1987_12 = {'days', 'wet weight'};

temp.tWw_HS1987_12  = C2K(data_T_RGR_HS1987(1,1)); units.temp.tWw_HS1987_12  = 'K'; label.temp.tWw_HS1987_12  = 'temperature'; 
light.tWw_HS1987_12 = 110 * 1e-6 * 3600; units.light.tWw_HS1987_12 = 'mol gamma m-2 s-1'; label.light.tWw_HS1987_12 = "light intensity"; 
nitrogen.tWw_HS1987_12 = 15 * 1e-6 ; units.nitrogen.tWw_HS1987_12 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_HS1987_12 = "N concentration"; 
CO2.tWw_HS1987_12 = 0.002 ; units.CO2.tWw_HS1987_12 = 'mol CO2 L-1'; label.CO2.tWw_HS1987_12 = "C concentration"; 
phosphorous.tWw_HS1987_12  =  15 * 1e-6; units.phosphorous.tWw_HS1987_12 = 'mol P043- L-1 '; label.phosphorous.tWw_HS1987_12 = "P concentration"; 
photoPeriod.tWw_HS1987_12   = 1; units.photoPeriod.tWw_HS1987_12 = '-'; label.photoPeriod.tWw_HS1987_12 = "Photoperiod"; 

bibkey.tWw_HS1987_12  = 'Hanisak-Samuel_1987'; 

% 18°C
Ww0.tWw_HS1987_18 = 26; units.Ww0.tWw_HS1987_18 = 'g Ww'; label.Ww0.tWw_HS1987_18 = 'Initial wet weight';
time.tWw_HS1987_18 = 21; units.time.tWw_HS1987_18 = 'd'; label.time.tWw_HS1987_18 = 'days';
WWf_HS1987_18 = 2^(data_T_RGR_HS1987(2,2) * time.tWw_HS1987_18) * Ww0.tWw_HS1987_18;
data.tWw_HS1987_18 = [0 Ww0.tWw_HS1987_18; time.tWw_HS1987_18 WWf_HS1987_18];
units.tWw_HS1987_18 = {'d', 'g Ww'}; label.tWw_HS1987_18 = {'days', 'wet weight'};

temp.tWw_HS1987_18  = C2K(data_T_RGR_HS1987(2,1)); units.temp.tWw_HS1987_18  = 'K'; label.temp.tWw_HS1987_18  = 'temperature'; 
light.tWw_HS1987_18 = 110 * 1e-6 * 3600; units.light.tWw_HS1987_18 = 'mol gamma m-2 s-1'; label.light.tWw_HS1987_18 = "light intensity"; 
nitrogen.tWw_HS1987_18 = 15 * 1e-6 ; units.nitrogen.tWw_HS1987_18 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_HS1987_18 = "N concentration"; 
CO2.tWw_HS1987_18 = 0.002 ; units.CO2.tWw_HS1987_18 = 'mol CO2 L-1'; label.CO2.tWw_HS1987_18 = "C concentration"; 
phosphorous.tWw_HS1987_18 = 15 * 1e-6; units.phosphorous.tWw_HS1987_18 = 'mol P043- L-1 '; label.phosphorous.tWw_HS1987_18 = "P concentration"; 
photoPeriod.tWw_HS1987_18 = 1; units.photoPeriod.tWw_HS1987_18 = '-'; label.photoPeriod.tWw_HS1987_18 = "Photoperiod"; 
bibkey.tWw_HS1987_18  = 'Hanisak-Samuel_1987';


% 24°C
Ww0.tWw_HS1987_24 = 26; units.Ww0.tWw_HS1987_24 = 'g Ww'; label.Ww0.tWw_HS1987_24 = 'Initial wet weight';
time.tWw_HS1987_24 = 21; units.time.tWw_HS1987_24 = 'd'; label.time.tWw_HS1987_24 = 'days';
WWf_HS1987_24 = 2^(data_T_RGR_HS1987(3,2) * time.tWw_HS1987_24) * Ww0.tWw_HS1987_24;
data.tWw_HS1987_24 = [0 Ww0.tWw_HS1987_24; time.tWw_HS1987_24 WWf_HS1987_24];
units.tWw_HS1987_24 = {'d', 'g Ww'}; label.tWw_HS1987_24 = {'days', 'wet weight'};


temp.tWw_HS1987_24  = C2K(data_T_RGR_HS1987(3,1)); units.temp.tWw_HS1987_24  = 'K'; label.temp.tWw_HS1987_24  = 'temperature'; 
light.tWw_HS1987_24 = 110 * 1e-6 * 3600; units.light.tWw_HS1987_24 = 'mol gamma m-2 s-1'; label.light.tWw_HS1987_24 = "light intensity"; 
nitrogen.tWw_HS1987_24 = 15 * 1e-6 ; units.nitrogen.tWw_HS1987_24 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_HS1987_24 = "N concentration"; 
CO2.tWw_HS1987_24 = 0.002 ; units.CO2.tWw_HS1987_24 = 'mol CO2 L-1'; label.CO2.tWw_HS1987_24 = "C concentration"; 
phosphorous.tWw_HS1987_24 = 15 * 1e-6; units.phosphorous.tWw_HS1987_24 = 'mol P043- L-1 '; label.phosphorous.tWw_HS1987_24 = "P concentration"; 
photoPeriod.tWw_HS1987_24 = 1; units.photoPeriod.tWw_HS1987_24 = '-'; label.photoPeriod.tWw_HS1987_24 = "Photoperiod"; 
bibkey.tWw_HS1987_24  = 'Hanisak-Samuel_1987';

% 30°C
Ww0.tWw_HS1987_30 = 26; units.Ww0.tWw_HS1987_30 = 'g Ww'; label.Ww0.tWw_HS1987_30 = 'Initial wet weight';
time.tWw_HS1987_30 = 21; units.time.tWw_HS1987_30 = 'd'; label.time.tWw_HS1987_30 = 'days';
WWf_HS1987_30 = 2^(data_T_RGR_HS1987(4,2) * time.tWw_HS1987_30) * Ww0.tWw_HS1987_30;
data.tWw_HS1987_30 = [0 Ww0.tWw_HS1987_30; time.tWw_HS1987_30 WWf_HS1987_30];
units.tWw_HS1987_30 = {'d', 'g Ww'}; label.tWw_HS1987_30 = {'days', 'wet weight'};

temp.tWw_HS1987_30  = C2K(data_T_RGR_HS1987(4,1)); units.temp.tWw_HS1987_30  = 'K'; label.temp.tWw_HS1987_30  = 'temperature'; 
light.tWw_HS1987_30 = 110 * 1e-6 * 3600; units.light.tWw_HS1987_30 = 'mol gamma m-2 s-1'; label.light.tWw_HS1987_30 = "light intensity"; 
nitrogen.tWw_HS1987_30 = 15 * 1e-6 ; units.nitrogen.tWw_HS1987_30 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_HS1987_30 = "N concentration"; 
CO2.tWw_HS1987_30 = 0.002 ; units.CO2.tWw_HS1987_30 = 'mol CO2 L-1'; label.CO2.tWw_HS1987_30 = "C concentration"; 
phosphorous.tWw_HS1987_30 = 15 * 1e-6; units.phosphorous.tWw_HS1987_30 = 'mol P043- L-1 '; label.phosphorous.tWw_HS1987_30 = "P concentration"; 
photoPeriod.tWw_HS1987_30 = 1; units.photoPeriod.tWw_HS1987_30 = '-'; label.photoPeriod.tWw_HS1987_30 = "Photoperiod"; 
bibkey.tWw_HS1987_30  = 'Hanisak-Samuel_1987';


% MG2023a
%time (days), total wet weight 
data_t_Ww_28 = [ ... % d, g Ww
                 0 12;
                 5 13.985915288639502;
                10 16.633802598053865;
                15 16.59154894752555;
                20 14.816901041131903];

%time (daysCNratio ), total wet weight (g Ww)
data_t_Ww_31 = [ ... % d, g Ww
                 0 12;
                 5 14.999999806578728;
                10 18.788732198674722;
                15 17.253520871589778;
                20 16.774647712972374 ]; 



data.tWw_MG2023a_28 = data_t_Ww_28(1:3, [1 2]); %Select 
units.tWw_MG2023a_28 = {'d', 'g Ww'}; label.tWw_MG2023a_28 = {'days' , 'wet weight'}; 
Ww0.tWw_MG2023a_28 = data_t_Ww_28(1,2); units.Ww0.tWw_MG2023a_28 = 'g Ww'; label.Ww0.tWw_MG2023a_31 = 'Initial wet weight'; 
time.tWw_MG2023a_28 = 10 ;  units.time.tWw_MG2023a_28 = 'd'; label.time.tWw_MG2023a_28 = 'days';


temp.tWw_MG2023a_28 = C2K(28); units.temp.tWw_MG2023a_28 = 'K'; label.temp.tWw_MG2023a_28 = 'temperature';
light.tWw_MG2023a_28 = mean([435 581]) * 1e-6 * 3600; units.light.tWw_MG2023a_28 = 'mol gamma m-2 s-1'; label.light.tWw_MG2023a_28 = "light intensity"; 
nitrogen.tWw_MG2023a_28 = (0.06 + 3.2) * 1e-6 ; units.nitrogen.tWw_MG2023a_28 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_MG2023a_28 = "N concentration"; 
CO2.tWw_MG2023a_28 = 0.002 ; units.CO2.tWw_MG2023a_28 = 'mol CO2 L-1'; label.CO2.tWw_MG2023a_28 = "C concentration"; 
phosphorous.tWw_MG2023a_28  =  0.5 * 1e-6; units.phosphorous.tWw_MG2023a_28 = 'mol P043- L-1 '; label.phosphorous.tWw_MG2023a_28 = "P concentration"; 
photoPeriod.tWw_MG2023a_28   = 1; units.photoPeriod.tWw_MG2023a_28 = '-'; label.photoPeriod.tWw_MG2023a_28 = "Photoperiod"; 

bibkey.tWw_MG2023a_28 = 'magana-gallegos_growth_2023'; 


data.tWw_MG2023a_31 = data_t_Ww_31(1:3, [1 2]); 
units.tWw_MG2023a_31 = {'d', 'g Ww'}; label.tWw_MG2023a_31 = {'days' , 'wet weight'}; 
Ww0.tWw_MG2023a_31 = data_t_Ww_31(1,2); units.Ww0.tWw_MG2023a_31 = 'g Ww'; label.Ww0.tWw_MG2023a_31 = 'Initial wet weight'; 
time.tWw_MG2023a_31 = 10 ;  units.time.tWw_MG2023a_31 = 'd'; label.time.tWw_MG2023a_31 = 'days';


temp.tWw_MG2023a_31 = C2K(31); units.temp.tWw_MG2023a_31 = 'K'; label.temp.tWw_MG2023a_31 = 'temperature';
light.tWw_MG2023a_31 = mean([435 581]) * 1e-6 * 3600; units.light.tWw_MG2023a_31 = 'mol gamma m-2 s-1'; label.light.tWw_MG2023a_31 = "light intensity"; 
nitrogen.tWw_MG2023a_31 = (0.06 + 3.2) * 1e-6 ; units.nitrogen.tWw_MG2023a_31 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_MG2023a_31 = "N concentration"; 
CO2.tWw_MG2023a_31 = 0.002 ; units.CO2.tWw_MG2023a_31 = 'mol CO2 L-1'; label.CO2.tWw_MG2023a_31 = "C concentration"; 
phosphorous.tWw_MG2023a_31  =  0.5 * 1e-6; units.phosphorous.tWw_MG2023a_31 = 'mol P043- L-1 '; label.phosphorous.tWw_MG2023a_31 = "P concentration"; 
photoPeriod.tWw_MG2023a_31   = 1; units.photoPeriod.tWw_MG2023a_31 = '-'; label.photoPeriod.tWw_MG2023a_31 = "Photoperiod"; 

bibkey.tWw_MG2023a_31 = 'magana-gallegos_growth_2023'; 




%MG2023b
% time (days), wet weight from growth rate (doubling d-1)
data_MG2023b_T_RGR = [  ... %°C ,   doubling d-1
                        22	0.078;
                        25  0.077;
                        28	0.095;
                        31  0.058];

Ww0.tWw_MG2023b_22 = 6; units.Ww0.tWw_MG2023b_22   = 'g Ww'; label.Ww0.tWw_MG2023b_22   = 'Initial wet weight';
time.tWw_MG2023b_22  = 5;  units.time.tWw_MG2023b_22  = 'd'; label.time.tWw_MG2023b_22   = 'days';
WWfMG2023b_22 = 2.^(data_MG2023b_T_RGR(1,2) * time.tWw_MG2023b_22 ) * Ww0.tWw_MG2023b_22 ;
data.tWw_MG2023b_22 = [0 Ww0.tWw_MG2023b_22 ; time.tWw_MG2023b_22   WWfMG2023b_22 ];
units.tWw_MG2023b_22 = {'d', 'g Ww'}; label.tWw_MG2023b_22 = {'days', 'wet weight'};

DGR_22 = ((WWfMG2023b_22 / Ww0.tWw_MG2023b_22)^(1/time.tWw_MG2023b_22) - 1) * 100; 
r_GR_22 = (1/time.tWw_MG2023b_22)* log(WWfMG2023b_22 / Ww0.tWw_MG2023b_22); 


temp.tWw_MG2023b_22  = C2K(data_MG2023b_T_RGR(1,1)); units.temp.tWw_MG2023b_22  = 'K'; label.temp.tWw_MG2023b_22  = 'temperature'; 
light.tWw_MG2023b_22 = mean([435 581]) * 1e-6 * 3600; units.light.tWw_MG2023b_22 = 'mol gamma m-2 s-1'; label.light.tWw_MG2023b_22 = "light intensity"; 
nitrogen.tWw_MG2023b_22 = (mean([0.06 0.05]) + mean([3.2 3.49])) * 1e-6 ; units.nitrogen.tWw_MG2023b_22 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_MG2023b_22 = "N concentration"; 
CO2.tWw_MG2023b_22 = 0.002 ; units.CO2.tWw_MG2023b_22 = 'mol CO2 L-1'; label.CO2.tWw_MG2023b_22 = "C concentration"; 
phosphorous.tWw_MG2023b_22  =  mean([0.57 0.5]) * 1e-6; units.phosphorous.tWw_MG2023b_22 = 'mol P043- L-1 '; label.phosphorous.tWw_MG2023b_22 = "P concentration"; 
photoPeriod.tWw_MG2023b_22   = 1; units.photoPeriod.tWw_MG2023b_22 = '-'; label.photoPeriod.tWw_MG2023b_22 = "Photoperiod"; 

bibkey.tWw_MG2023b_22  = 'magana-gallegos_temperature_2023'; 




Ww0.tWw_MG2023b_25 = 6; units.Ww0.tWw_MG2023b_25   = 'g Ww'; label.Ww0.tWw_MG2023b_25   = 'Initial wet weight';
time.tWw_MG2023b_25  = 5;  units.time.tWw_MG2023b_25  = 'd'; label.time.tWw_MG2023b_25   = 'days';
WWfMG2023b_25 = 2.^(data_MG2023b_T_RGR(2,2) * time.tWw_MG2023b_25 ) * Ww0.tWw_MG2023b_25;
data.tWw_MG2023b_25 = [0 Ww0.tWw_MG2023b_25 ; time.tWw_MG2023b_25   WWfMG2023b_25 ];
units.tWw_MG2023b_25 = {'d', 'g Ww'}; label.tWw_MG2023b_25 = {'days', 'wet weight'};

DGR_25 = ((WWfMG2023b_25 / Ww0.tWw_MG2023b_25)^(1/time.tWw_MG2023b_25) - 1) * 100; 
r_GR_25 = (1/time.tWw_MG2023b_25)* log(WWfMG2023b_25 / Ww0.tWw_MG2023b_25); 


temp.tWw_MG2023b_25  = C2K(data_MG2023b_T_RGR(2,1)); units.temp.tWw_MG2023b_25  = 'K'; label.temp.tWw_MG2023b_25  = 'temperature'; 
light.tWw_MG2023b_25 = mean([435 581]) * 1e-6 * 3600; units.light.tWw_MG2023b_25 = 'mol gamma m-2 s-1'; label.light.tWw_MG2023b_25 = "light intensity"; 
nitrogen.tWw_MG2023b_25 = (mean([0.06 0.05]) + mean([3.2 3.49])) * 1e-6 ; units.nitrogen.tWw_MG2023b_25 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_MG2023b_25 = "N concentration"; 
CO2.tWw_MG2023b_25 = 0.002 ; units.CO2.tWw_MG2023b_25 = 'mol CO2 L-1'; label.CO2.tWw_MG2023b_25 = "C concentration"; 
phosphorous.tWw_MG2023b_25  =  mean([0.57 0.5]) * 1e-6; units.phosphorous.tWw_MG2023b_25 = 'mol P043- L-1 '; label.phosphorous.tWw_MG2023b_25 = "P concentration"; 
photoPeriod.tWw_MG2023b_25   = 1; units.photoPeriod.tWw_MG2023b_25 = '-'; label.photoPeriod.tWw_MG2023b_25 = "Photoperiod"; 


bibkey.tWw_MG2023b_25  = 'magana-gallegos_temperature_2023'; 



Ww0.tWw_MG2023b_28 = 6; units.Ww0.tWw_MG2023b_28   = 'g Ww'; label.Ww0.tWw_MG2023b_28   = 'Initial wet weight';
time.tWw_MG2023b_28  = 5;  units.time.tWw_MG2023b_28  = 'd'; label.time.tWw_MG2023b_28   = 'days';
WWfMG2023b_28 = 2.^(data_MG2023b_T_RGR(3,2) * time.tWw_MG2023b_28 ) * Ww0.tWw_MG2023b_28;
data.tWw_MG2023b_28 = [0 Ww0.tWw_MG2023b_28 ; time.tWw_MG2023b_28   WWfMG2023b_28 ];
units.tWw_MG2023b_28 = {'d', 'g Ww'}; label.tWw_MG2023b_28 = {'days', 'wet weight'};

DGR_28 = ((WWfMG2023b_28 / Ww0.tWw_MG2023b_28)^(1/time.tWw_MG2023b_28) - 1) * 100; 
r_GR_28 = (1/time.tWw_MG2023b_28)* log(WWfMG2023b_28 / Ww0.tWw_MG2023b_28); 

temp.tWw_MG2023b_28  = C2K(data_MG2023b_T_RGR(3,1)); units.temp.tWw_MG2023b_28  = 'K'; label.temp.tWw_MG2023b_28  = 'temperature'; 
light.tWw_MG2023b_28 = mean([435 581]) * 1e-6 * 3600; units.light.tWw_MG2023b_28 = 'mol gamma m-2 s-1'; label.light.tWw_MG2023b_28 = "light intensity"; 
nitrogen.tWw_MG2023b_28 = (mean([0.06 0.05]) + mean([3.2 3.49])) * 1e-6 ; units.nitrogen.tWw_MG2023b_28 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_MG2023b_28 = "N concentration"; 
CO2.tWw_MG2023b_28 = 0.002 ; units.CO2.tWw_MG2023b_28 = 'mol CO2 L-1'; label.CO2.tWw_MG2023b_28 = "C concentration"; 
phosphorous.tWw_MG2023b_28  = mean([0.57 0.5]) * 1e-6; units.phosphorous.tWw_MG2023b_28 = 'mol P043- L-1 '; label.phosphorous.tWw_MG2023b_28 = "P concentration"; 
photoPeriod.tWw_MG2023b_28   = 1; units.photoPeriod.tWw_MG2023b_28 = '-'; label.photoPeriod.tWw_MG2023b_28 = "Photoperiod"; 


bibkey.tWw_MG2023b_28  = 'magana-gallegos_temperature_2023'; 



Ww0.tWw_MG2023b_31 = 6; units.Ww0.tWw_MG2023b_31   = 'g Ww'; label.Ww0.tWw_MG2023b_31   = 'Initial wet weight';
time.tWw_MG2023b_31  = 5;  units.time.tWw_MG2023b_31  = 'd'; label.time.tWw_MG2023b_31   = 'days';
WWfMG2023b_31 = 2.^(data_MG2023b_T_RGR(4,2) * time.tWw_MG2023b_31) * Ww0.tWw_MG2023b_31 ;
data.tWw_MG2023b_31 = [0 Ww0.tWw_MG2023b_31 ; time.tWw_MG2023b_31  WWfMG2023b_31];
units.tWw_MG2023b_31 = {'d', 'g Ww'}; label.tWw_MG2023b_31 = {'days', 'wet weight'};

DGR_31 = ((WWfMG2023b_31 / Ww0.tWw_MG2023b_31)^(1/time.tWw_MG2023b_31) - 1) * 100; 
r_GR_31 = (1/time.tWw_MG2023b_31)* log(WWfMG2023b_31 / Ww0.tWw_MG2023b_31); 

temp.tWw_MG2023b_31  = C2K(data_MG2023b_T_RGR(4,1)); units.temp.tWw_MG2023b_31  = 'K'; label.temp.tWw_MG2023b_31  = 'temperature'; 
light.tWw_MG2023b_31 = mean([435 581]) * 1e-6 * 3600; units.light.tWw_MG2023b_31 = 'mol gamma m-2 s-1'; label.light.tWw_MG2023b_31 = "light intensity"; 
nitrogen.tWw_MG2023b_31 = (mean([0.06 0.05]) + mean([3.2 3.49])) * 1e-6 ; units.nitrogen.tWw_MG2023b_31 = 'mol NO3 + NH4 L-1'; label.nitrogen.tWw_MG2023b_31 = "N concentration"; 
CO2.tWw_MG2023b_31 = 0.002 ; units.CO2.tWw_MG2023b_31 = 'mol CO2 L-1'; label.CO2.tWw_MG2023b_31 = "C concentration"; 
phosphorous.tWw_MG2023b_31  =  mean([0.57 0.5]) * 1e-6; units.phosphorous.tWw_MG2023b_31 = 'mol P043- L-1 '; label.phosphorous.tWw_MG2023b_31 = "P concentration"; 
photoPeriod.tWw_MG2023b_31   = 1; units.photoPeriod.tWw_MG2023b_31 = '-'; label.photoPeriod.tWw_MG2023b_31 = "Photoperiod"; 


bibkey.tWw_MG2023b_31  = 'magana-gallegos_temperature_2023'; 






% Irradiance - O2 photosynthesis rate 

data_I_J_O2_VE2023_HL = [ ... % Light intensity [micro mol quanta m-2 s-1], gross photosynthesis rate [micro mol O2 cm-2 h-1] %Data obtained from Roman               
                       0	0
                    12.5	0.322841033
                      25	0.388713675
                      44	0.561393453
                      65	0.725846695
                      85	0.832121022
                      99	0.967540006
                      150	1.081253676
                      185	1.171363805
                      346	1.329566572
                      648	1.408211531
                      994	1.40457139
                      1228	1.389586446
                      1500	1.466176878];

dwratio = 0.178;      
smcoeff = 0.1727;    

data.I_J_O2_VE2023_HL(:,1) =  data_I_J_O2_VE2023_HL(:,1); 
data.I_J_O2_VE2023_HL(:,2) =  data_I_J_O2_VE2023_HL(:,2); 
% data.I_J_O2_VE2023_HL(:,2) =  data_I_J_O2_VE2023_HL(:,2) / smcoeff/ dwratio ; %micro mol O2 gdw-1 h-1 For comparaison with Venolia 
units.I_J_O2_VE2023_HL = {'micro mol quanta m-2 s-1', 'micro mol O2 cm-2 h-1'};
% units.I_J_O2_VE2023_HL = {'micro mol quanta m-2 s-1', 'micro mol O2 gdW-1
% h-1'};w
label.I_J_O2_VE2023_HL = {'Light intensity', 'Gross photosynthesis rate'}; 

surface.I_J_O2_VE2023_HL = 1.5 ; units.surface.I_J_O2_VE2023_HL = 'cm^2'; label.surface.I_J_O2_VE2023_HL = 'Initial surface'; 
temp.I_J_O2_VE2023_HL = C2K(23); units.temp.I_J_O2_VE2023_HL = 'K'; label.temp.I_J_O2_VE2023_HL = 'Temperature'; 
photoPeriod.I_J_O2_VE2023_HL   = 1; units.photoPeriod.I_J_O2_VE2023_HL = '-'; label.photoPeriod.I_J_O2_VE2023_HL = "Photoperiod"; 


bibkey.I_J_O2_VE2023_HL =  'vasquezelizondo_growth_2023'; 
comment.I_J_O2_VE2023_HL = 'SF3 adapted to high light conditions before testing irradiance effect (HL corresponds to = 657 +- 11 micro mol quanta m-2 s-1 at 12:12 period ';




% data_I_J_O2_VE2023_LL = [ ...  % Light intensity [micro mol quanta m-2 s-1], gross photosynthesis rate [micro mol O2 cm-2 h-1]  %Data obtained from Roman
%                         0       0
%                         25      0.49317336
%                         44      0.703456499
%                         65      0.813412543
%                         85      0.848097986
%                         99      0.880773596
%                         150     0.91267179
%                         185     0.916707401
%                         346     0.897195079
%                         648     0.842498429
%                         944     0.8337302
%                         1228    0.779897532
%                         1504    0.759838253];
% 
% 
% data.I_J_O2_VE2023_LL(:,1) =  data_I_J_O2_VE2023_LL(:,1) ; % 
% % data.I_J_O2_VE2023_LL(:,2) =  data_I_J_O2_VE2023_LL(:,2)  ; % 
% data.I_J_O2_VE2023_LL(:,2) =  data_I_J_O2_VE2023_LL(:,2) / smcoeff/ dwratio ; % 
% % units.I_J_O2_VE2023_LL = {'micro mol quanta m-2 s-1', 'micro mol O2 cm-2 h-1'};
% units.I_J_O2_VE2023_LL  = {'micro mol quanta m-2 s-1', 'micro mol O2 gdW-1 h-1'};
% label.I_J_O2_VE2023_LL= {'Light intensity', 'Gross photosynthesis rate'}; 
% 
% surface.I_J_O2_VE2023_LL = 1.5 ; units.surface.I_J_O2_VE2023_LL = 'cm^2'; label.surface.I_J_O2_VE2023_LL = 'Initial surface'; 
% temp.I_J_O2_VE2023_LL = C2K(23); units.temp.I_J_O2_VE2023_LL = 'K'; label.temp.I_J_O2_VE2023_LL = 'Temperature'; 
% photoPeriod.I_J_O2_VE2023_LL   = 1; units.photoPeriod.I_J_O2_VE2023_LL = '-'; label.photoPeriod.I_J_O2_VE2023_LL = "Photoperiod"; 
% 
% bibkey.I_J_O2_VE2023_LL =  'vasquezelizondo_growth_2023'; 
% comment.I_J_O2_VE2023_LL = 'SF3 adapted to low light conditions before testing irradiance effect (LL corresponds to = 105 +-4 micro mol quant m-1 s-1 at 12:12 period';
% 
% 
% 

%% New datasets 

if contains(pwd,'Identifiability')
    load("dataExperiment.mat"); 
    load("dataAvailable.mat");
    data = dataAvailable; 
    
    nm = fieldnames(dataExperiment); 
    nst = length(nm); 
    
    for i =  1:nst
        nm_exp = string(nm(i));
        if contains(nm_exp,'TC_Ww')
            data.(nm_exp) = [dataExperiment.(nm_exp).timePoints',  dataExperiment.(nm_exp).Ww] ;
            data.(nm_exp) = [dataExperiment.(nm_exp).timePoints',  dataExperiment.(nm_exp).Ww] ;
    
            Ww0.(nm_exp) =  dataExperiment.(nm_exp).Ww0;  units.Ww0.(nm_exp)   = 'g Ww'; label.Ww0.(nm_exp)   = 'Initial wet weight';
            time.(nm_exp)  = dataExperiment.(nm_exp).time;  units.time.(nm_exp)  = 'd'; label.time.(nm_exp)   = 'days';
            units.(nm_exp) = {'d', 'g Ww'}; label.(nm_exp) = {'days', 'wet weight'};
    
    
            temp.(nm_exp)   = dataExperiment.(nm_exp).temp; units.temp.(nm_exp)   = 'K'; label.temp.(nm_exp)   = 'temperature'; 
            light.(nm_exp)  = dataExperiment.(nm_exp).light; units.light.(nm_exp)  = 'mol gamma m-2 s-1'; label.light.(nm_exp)  = "light intensity"; 
            nitrogen.(nm_exp)  =  dataExperiment.(nm_exp).nitrogen ; units.nitrogen.(nm_exp)  = 'mol NO3 + NH4 L-1'; label.nitrogen.(nm_exp)  = "N concentration"; 
            CO2.(nm_exp)  = dataExperiment.(nm_exp).CO2 ; units.CO2.(nm_exp)  = 'mol CO2 L-1'; label.CO2.(nm_exp) = "C concentration"; 
            phosphorous.(nm_exp)   = dataExperiment.(nm_exp).phosphorous; units.phosphorous.(nm_exp)  = 'mol P043- L-1 '; label.phosphorous.(nm_exp)  = "P concentration"; 
            photoPeriod.(nm_exp)   = dataExperiment.(nm_exp).photoPeriod; units.photoPeriod.(nm_exp)  = '- '; label.photoPeriod.(nm_exp)  = "Photoperiod"; 
    
            bibkey.(nm_exp)  = nm_exp; 
    
    
        elseif contains(nm_exp,'Starvation')
            % data.(nm_exp) = [dataExperiment.(nm_exp).timePoints',  dataExperiment.(nm_exp).Ww] ;
            data_Ww_Starvation = dataExperiment.(nm_exp).Ww; 
            data_CNratio_Starvation = dataExperiment.(nm_exp).CNtotalRatio; 

            data.(sprintf('tWw_%s', nm_exp))= [dataExperiment.(nm_exp).timePoints', data_Ww_Starvation];
            
            data.(sprintf('CNratio_%s', nm_exp)) = [dataExperiment.(nm_exp).timePoints', data_CNratio_Starvation];

            Ww0.(sprintf('tWw_%s', nm_exp)) =  dataExperiment.(nm_exp).Ww0;  units.Ww0.(sprintf('tWw_%s', nm_exp))   = 'g Ww'; label.Ww0.(sprintf('tWw_%s', nm_exp))  = 'Initial wet weight';
            time.(nm_exp)  = dataExperiment.(nm_exp).time;  units.time.(nm_exp)  = 'd'; label.time.(nm_exp)   = 'days';
            
            units.(sprintf('tWw_%s', nm_exp)) = {'d', 'g Ww'}; label.(sprintf('tWw_%s', nm_exp)) = {'days', 'wet weight'};
            units.(sprintf('CNratio_%s', nm_exp)) = {'d', 'mol C mol N'}; label.(sprintf('CNratio_%s', nm_exp)) = {'days','CN ratio'};
    
    
            temp.(nm_exp)   = dataExperiment.(nm_exp).temp; units.temp.(nm_exp)   = 'K'; label.temp.(nm_exp)   = 'temperature'; 
            light.(nm_exp)  = dataExperiment.(nm_exp).light; units.light.(nm_exp)  = 'mol gamma m-2 s-1'; label.light.(nm_exp)  = "light intensity"; 
            nitrogen.(nm_exp)  =  dataExperiment.(nm_exp).nitrogen ; units.nitrogen.(nm_exp)  = 'mol NO3 + NH4 L-1'; label.nitrogen.(nm_exp)  = "N concentration"; 
            CO2.(nm_exp)  = dataExperiment.(nm_exp).CO2 ; units.CO2.(nm_exp)  = 'mol CO2 L-1'; label.CO2.(nm_exp) = "C concentration"; 
            phosphorous.(nm_exp)   = dataExperiment.(nm_exp).phosphorous; units.phosphorous.(nm_exp)  = 'mol P043- L-1 '; label.phosphorous.(nm_exp)  = "P concentration"; 
            photoPeriod.(nm_exp)   = dataExperiment.(nm_exp).photoPeriod; units.photoPeriod.(nm_exp)  = '- '; label.photoPeriod.(nm_exp)  = "Photoperiod"; 

        elseif contains(nm_exp,'N_uptake')
            data.(nm_exp) = [(dataExperiment.(nm_exp).NitrateConcentration * 1e-6)',  (dataExperiment.(nm_exp).NitrateAssimilation * 1e-6)'] ;
            Ww0.(nm_exp) =  dataExperiment.(nm_exp).Ww0;  units.Ww0.(nm_exp)   = 'g Ww'; label.Ww0.(nm_exp)   = 'Initial wet weight';
            units.(nm_exp) = {'mol NO3 L-1', 'mol N g dW-1 hr-1'}; label.(nm_exp) = {'Nitrate concentration', 'Nitrogen uptake rate'}; 
            temp.(nm_exp)   = C2K(dataExperiment.(nm_exp).temp); units.temp.(nm_exp)   = 'K'; label.temp.(nm_exp)   = 'temperature'; 
            
        end
    end
end


%%  set weights for all real data
weights = setweights(data, []); %% 

weights.tWw_MG2023a_28 =  weights.tWw_MG2023a_28 * 0.25; %Reduce weight as is the opposite from MG
weights.tWw_MG2023a_31 =  weights.tWw_MG2023a_31 * 0.25;

weights.tWw_MG2023b_22 = weights.tWw_MG2023b_22 * 10; 
weights.tWw_MG2023b_25 = weights.tWw_MG2023b_25 * 10; 
weights.tWw_MG2023b_28 = weights.tWw_MG2023b_28 * 10; 
weights.tWw_MG2023b_31 = weights.tWw_MG2023b_31 * 10; 


weights.I_J_O2_VE2023_HL=  weights.I_J_O2_VE2023_HL * 0.25;
% weights.I_J_O2_VE2023_LL =  weights.I_J_O2_VE2023_LL * 1;


weights.tWw_HS1987_12 = weights.tWw_HS1987_12 * 0; 
weights.tWw_HS1987_18 = weights.tWw_HS1987_18 * 0; 
weights.tWw_HS1987_24 = weights.tWw_HS1987_24 * 0; 
weights.tWw_HS1987_30 = weights.tWw_HS1987_30 * 0; 





%% Setting up weights and auxdata if in Identifiability analysis

if contains(pwd,'Identifiability')
weights = setweights(data, []); %% 
    for i =  1:nst
        nm_exp = string(nm(i));
         if contains(nm_exp,'Starvation')
            weights.(sprintf('tWw_%s', nm_exp))= weights.(sprintf('tWw_%s', nm_exp));
            weights.(sprintf('CNratio_%s', nm_exp)) = weights.(sprintf('CNratio_%s', nm_exp));
            auxData.Ww0.(sprintf('tWw_%s', nm_exp)) =  Ww0.(sprintf('tWw_%s', nm_exp)) ;
            if isfield(time,(nm_exp))
                auxData.time.(sprintf('tWw_%s', nm_exp)) =  time.(nm_exp) ;
             end
         else
             weights.(nm_exp) =  weights.(nm_exp);
             auxData.Ww0.(nm_exp) =  Ww0.(nm_exp) ;
             if isfield(time,(nm_exp))
                auxData.time.(nm_exp) =  time.(nm_exp) ;
             end
         end
    end
end


%% pack auxData and txtData for output
auxData.temp = temp;
auxData.light = light;
auxData.nitrogen = nitrogen;
auxData.CO2 = CO2;
auxData.phosphorous = phosphorous;
auxData.photoPeriod = photoPeriod;


auxData.Ww0.HT2023_CNratio = Ww0.HT2023_CNratio;


auxData.Ww0.tWw_MG2023a_28 = Ww0.tWw_MG2023a_28;
auxData.Ww0.tWw_MG2023a_31 = Ww0.tWw_MG2023a_31;

auxData.Ww0.tWw_MG2023b_22 = Ww0.tWw_MG2023b_22;
auxData.Ww0.tWw_MG2023b_25 = Ww0.tWw_MG2023b_25;
auxData.Ww0.tWw_MG2023b_28 = Ww0.tWw_MG2023b_28;
auxData.Ww0.tWw_MG2023b_31 = Ww0.tWw_MG2023b_31;

auxData.Ww0.tWw_HS1987_12 = Ww0.tWw_HS1987_12;
auxData.Ww0.tWw_HS1987_18 = Ww0.tWw_HS1987_18;
auxData.Ww0.tWw_HS1987_24 = Ww0.tWw_HS1987_24;
auxData.Ww0.tWw_HS1987_30 = Ww0.tWw_HS1987_30;




auxData.surface.I_J_O2_VE2023_HL = surface.I_J_O2_VE2023_HL; 
% auxData.surface.I_J_O2_VE2023_LL = surface.I_J_O2_VE2023_LL ; 



auxData.time.HT2023_CNratio = time.HT2023_CNratio;

auxData.time.tWw_MG2023a_28 = time.tWw_MG2023a_28;
auxData.time.tWw_MG2023a_31 = time.tWw_MG2023a_31;

auxData.time.tWw_MG2023b_22 = time.tWw_MG2023b_22;
auxData.time.tWw_MG2023b_25 = time.tWw_MG2023b_25;
auxData.time.tWw_MG2023b_28 = time.tWw_MG2023b_28;
auxData.time.tWw_MG2023b_31 = time.tWw_MG2023b_31;


auxData.time.tWw_HS1987_12 = time.tWw_HS1987_12;
auxData.time.tWw_HS1987_18 = time.tWw_HS1987_18;
auxData.time.tWw_HS1987_24 = time.tWw_HS1987_24;
auxData.time.tWw_HS1987_30 = time.tWw_HS1987_30;



txtData.units = units; 
txtData.label = label;
txtData.bibkey = bibkey;


%% Group plots 
set1 = {'tWw_MG2023a_28','tWw_MG2023a_31'}; subtitle1 = {'Data of wet weight MG2023a'};
set2 = {'tWw_MG2023b_22','tWw_MG2023b_25', 'tWw_MG2023b_28', 'tWw_MG2023b_31'}; subtitle2 = {'Data of wet weight MG2023b'};
set3 = {'tWw_HS1987_12','tWw_HS1987_18','tWw_HS1987_24','tWw_HS1987_30'}; subtitle3 = {'Data of HS1987'};  


metaData.grp.sets = {set1,set2,set3};
metaData.grp.subtitle = {subtitle1, subtitle2, subtitle3};


 
if contains(pwd,'Identifiability')
    set4 = {}; 
    set5 = {}; 
    set6 = {}; 

    for i =  1:nst
        nm_exp = char(nm(i));
        if contains(nm_exp, 'TC_Ww_') && ~contains(nm_exp, 'twoDays')
            set4{end+1} = nm_exp;
            subtitle4 = {'Data of TC_Ww'}; 
        elseif contains(nm_exp, 'N_uptake')
            set5{end+1} = nm_exp;
            subtitle5 = {'Data of N_uptake'}; 
        elseif contains(nm_exp, 'TC_Ww_twoDays_')
            set6{end+1} = nm_exp; 
            subtitle6 = {'Data of TC_Ww_twoDays'}; 
        end
    end
metaData.grp.sets = {set1,set2,set3, set4, set5, set6};
metaData.grp.subtitle = {subtitle1, subtitle2,subtitle3,subtitle4, subtitle5, subtitle6};

end

