function [data, auxData, metaData, txtData, weights] = mydata_Sargassum_fluitans_Arrhenius

%% set metaData
metaData.phylum     = 'Ochrophyta'; 
metaData.class      = 'Phaeophycaeae'; 
metaData.order      = 'Fucales'; 
metaData.family     = 'Sargassaceae';
metaData.species    = 'Sargassum_fluitans_Arrhenius'; 
metaData.species_en = 'SFF'; 


metaData.T_typical  = C2K(25); % K, 

metaData.data_1     = {'T-RGR'; 't-Ww'; 'I-J_O2'}; 

%% set data
% uni-variate data
%% set data
%% Import from tempcorr factor
load("data_TC_GR.mat"); 

%% Changing values by %  
vectorData(1:4,2) =  vectorData(1:4,2) + (vectorData(1:4,2) * 0.2); % Modify value at 28°C to be higher
vectorData(8,2)=  vectorData(8,2) + (vectorData(8,2)  * 0.1); % Modify value at 28°C to be higher

data_vector_GR = vectorData; 

%% Creating data structure

data.TC_RG = [data_vector_GR(:,1) data_vector_GR(:,2) ];
units.TC_RG = {'T', '-'}; label.TC_RG = {'°C', 'ct'};

temp.TC_RG  = C2K(data_vector_GR(:,1)); units.temp.TC_RG = 'K'; label.temp.TC_RG  = 'temperature'; 


bibkey.TC_RG = 'TC_RG'; 


%%  set weights for all real data


weights = setweights(data, []); %% 
weights.TC_RG(8,:) = weights.TC_RG(8,:) * 1.5; 
%set pseudodata
data.psd.T_A = 3000;     units.psd.T_A = 'K';       label.psd.T_A = 'Arrhenius temperature';
data.psd.T_H = C2K(32);     units.psd.T_H = 'K';       label.psd.T_H = 'Arrhenius temperature';
data.psd.T_L = C2K(18);     units.psd.T_L = 'K';       label.psd.T_L = 'Arrhenius temperature';
data.psd.T_AH = 60000;     units.psd.T_AH = 'K';       label.psd.T_AH = 'Arrhenius temperature';
data.psd.T_AL = 40000;     units.psd.T_AL = 'K';       label.psd.T_AL = 'Arrhenius temperature';


weights.psd = setweights(data.psd, []);
weights.psd.T_A   = 0 * weights.psd.T_A;
weights.psd.T_H  = 0 * weights.psd.T_H;
weights.psd.T_L   = 0 * weights.psd.T_L;
weights.psd.T_AH   = 0 * weights.psd.T_AH;
weights.psd.T_AL   = 0 * weights.psd.T_AL;





%% pack auxData and txtData for output

auxData.temp = temp;



txtData.units = units; 
txtData.label = label;
txtData.bibkey = bibkey;
