%% -------------------- DATA on SARGASSUM COMPOSITION ---------------------%%
%% Dry to weight ratio
data_SAVEC_raw = readtable("Ratio MS et MF Sargassum_SAVE-C.xlsx"); 
data_SAVEC_raw.Station = string(data_SAVEC_raw.Station);
data_SAVEC_nostranding = data_SAVEC_raw(data_SAVEC_raw.Station ~= "Echouement", :);

data_SAVEC_nostranding_SFF = data_SAVEC_nostranding(data_SAVEC_nostranding.Morphotype == "S3", :);
mean_dwratio_SFF = mean(data_SAVEC_nostranding_SFF.MS_MF); 

data_SAVEC_nostranding_SNN = data_SAVEC_nostranding(data_SAVEC_nostranding.Morphotype == "S1", :);
mean_dwratio_SNN = mean(data_SAVEC_nostranding_SFF.MS_MF); 

data_SAVEC_nostranding_SNW = data_SAVEC_nostranding(data_SAVEC_nostranding.Morphotype == "S2", :);
mean_dwratio_SNW = mean(data_SAVEC_nostranding_SFF.MS_MF); 


%% Compound composition
% TON2021
data_TON2021_raw = readtable('CompositionSargassum - Tononetal_2021_ML.csv'); 
data_TON2021_compounds = unique(data_TON2021_raw.Compound_unit_);
mean_TON2021_compounds = struct();

for i = 1:length(data_TON2021_compounds)
    compoundName = data_TON2021_compounds{i,:}; 
    rows = strcmp(data_TON2021_raw.Compound_unit_, compoundName);
    values = data_TON2021_raw.Mean(rows);
    fieldName_compound = matlab.lang.makeValidName(compoundName);
    mean_TON2021_compounds.(fieldName_compound) = mean(values, 'omitNaN');
end

% MCH2022
data_MCH2022_raw = readtable('CompositionSargassum - Machadoetal2022_ML.csv'); 
data_MCH2022_compounds = unique(data_MCH2022_raw.BiochemicalContent);
mean_MCH2022_compounds = struct();

for i = 1:length(data_MCH2022_compounds)
    compoundName = data_MCH2022_compounds{i,:}; 
    rows = strcmp(data_MCH2022_raw.BiochemicalContent, compoundName);
    values = data_MCH2022_raw.x_DW(rows);
    fieldName_compound = matlab.lang.makeValidName(compoundName);
    mean_MCH2022_compounds.(fieldName_compound) = mean(values, 'omitNaN');
end


% MCH2024
data_MCH2024_raw = readtable('CompositionSargassum - Machadoetal_2024_ML.csv'); 

data_MCH2024_Moisture = data_MCH2024_raw.Moisture__BiomassDW_;
mean_MCH2024_moisture = mean(data_MCH2024_Moisture); 
data_MCH2024_Ash = data_MCH2024_raw.Ash__BiomassDW_;
mean_MCH2024_ash = mean(data_MCH2024_Ash); 
data_MCH2024_Proteins = data_MCH2024_raw.Proteins__BiomassDW_;
mean_MCH2024_Proteins = mean(data_MCH2024_Proteins); 
data_MCH2024_Phlorotannins = data_MCH2024_raw.Phlorotannins__BiomassDW_;
mean_MCH2024_Phlorotannins = mean(data_MCH2024_Phlorotannins); 
data_MCH2024_fucoxanthins = data_MCH2024_raw.fucoxanthin__BiomassDW;
mean_MCH2024_fucoxanthins = mean(data_MCH2024_fucoxanthins); 

mean_MCH2024_compounds = struct( ...
    'Moisture__BiomassDW_', mean_MCH2024_moisture, ...
    'Ash__BiomassDW_', mean_MCH2024_ash, ...
    'Proteins__BiomassDW_', mean_MCH2024_Proteins, ...
    'Phlorotannins__BiomassDW_', mean_MCH2024_Phlorotannins, ...
    'fucoxanthin__BiomassDW', mean_MCH2024_fucoxanthins);


%%
mean_TON2021_compounds = changeFieldNames(mean_TON2021_compounds);
mean_MCH2022_compounds = changeFieldNames(mean_MCH2022_compounds);
mean_MCH2024_compounds = changeFieldNames(mean_MCH2024_compounds);
datasets = {mean_TON2021_compounds, mean_MCH2022_compounds, mean_MCH2024_compounds};  % Example

allFields = {};
for i = 1:length(datasets)
    allFields = [allFields, fieldnames(datasets{i})'];
end

D = regexp(allFields,'[^_]+','once','match');
uniqueFields = unique(D);

% Calculate the mean values
meanValues = struct(); 

for i = 1:length(uniqueFields)
    field = uniqueFields{i};
    values = [];
    for j = 1:length(datasets)
        if isfield(datasets{j}, field)
            values(end+1) = datasets{j}.(field);
        end
    end

    if ~isempty(values)
        meanValues.(field) = mean(values);  % Mean across available datasets
    else
        meanValues.(field) = NaN;  % <-- Optional: fill with NaN if missing everywhere
    end
end



%% Plot compounds propotions per compound
meanVector_Compounds = struct2array(meanValues);

bar(1,meanVector_Compounds, 'stacked')
legend(uniqueFields)


%% Extract per organic inorganic


Wd_percentage_SFF = 100 * mean_dwratio_SFF; 
Ww_percentage_SFF =  100 - Wd_percentage_SFF; 

originalStruct = meanValues;

groupedCompositionStruct = struct();
groupedCompositionStruct.Ash = originalStruct.Ash;
groupedCompositionStruct.Moisture = originalStruct.Moisture;
groupedCompositionStruct.OrganicKnown = 0;
groupedCompositionStruct.OrganicNotKnown = 0;

structFieldNames  = fieldnames(meanValues);
InorganicMatterNames = {'Ash', 'Moisture'};

for i = 1:length(structFieldNames)
    field = structFieldNames{i};
    if ~ismember(field, InorganicMatterNames)
        groupedCompositionStruct.OrganicKnown = groupedCompositionStruct.OrganicKnown + originalStruct.(field);   
    end
    
end


groupedCompositionStruct.OrganicNotKnown = 100 - sum(struct2array(groupedCompositionStruct));
%% Plot compounds
figure
t = tiledlayout(1,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
ah = bar(categorical({'Sargassum spp.'}), [Wd_percentage_SFF Ww_percentage_SFF ], 'stacked'); 
set(ah, 'FaceColor', 'Flat')
ah(1).CData = min([ 0.9023    0.5613    0.1044] ,1) ;  
ah(2).CData = min([0.5184    0.7584    0.7230] ,1) ;  

lg = legend('Dry weight', 'Water', 'Orientation', 'vertical');
    
lg.Layout.Tile = 'north';
ylabel('%Total biomass')


nexttile
bh = bar(categorical({'General compound Sargassum'}),struct2array(groupedCompositionStruct), 'stacked'); 
set(bh, 'FaceColor', 'Flat')

bh(1).CData = min([0.2196    0.0078    0.2314] + 0.4,1) ;  
bh(2).CData = min([0.0081    0.4822    0.6031] + 0.4,1) ;  
bh(3).CData = min([0.9391    0.3499    0.0103] + 0.4,1) ;  
bh(4).CData = min([0.6377    0.0505    0.1249] + 0.4,1) ;   

ylabel('%biomass gDw')

lg = legend(unique(fieldnames(groupedCompositionStruct)), 'Orientation', 'vertical');
    
lg.Layout.Tile = 'north';
%% 
hAx = axes; 
hAx.FontSize = 30;
%% Exp compounds still to do --> order ash moisture and the rest
% nexttile
% ch = bar(categorical({'Compound Sargassum'}),meanVector_Compounds, 'stacked'); 
% set(ch, 'FaceColor', 'Flat')
% 
% legend(uniqueFields)
% % 
% 
% ch(1).CData = min([0.2196    0.0078    0.2314] + 0.4,1) ;  
% ch(2).CData = min([0.0081    0.4822    0.6031] + 0.4,1) ; 
% for i = 3:length(ch)
%     ch(i).CData = min(([0.9391    0.3499    0.0103] + (i/10)),1) ;  
% end
% 
% % 

%% Elemental composition
% Mean Sargassum
mean_MG2023_elementalComposition = struct( ...
    'Carbon', 27.448633014, ...
    'Nitrogen', 1.221088435, ...
    'Phosphorous', 0.08065068493);

% HTT2024
data_HTT2024_raw = readtable('CompositionSargassum - Lapointeetal_2021_ML.csv'); 
data_HTT2024_elements = unique(data_HTT2024_raw.Var2);
mean_HTT2024_elements = struct();


for i = 1:length(data_HTT2024_elements)
    elementName = data_HTT2024_elements{i,:}; 
    rows = strcmp(data_HTT2024_raw.Var2, elementName);
    values = data_HTT2024_raw.SFF(rows);
    fieldName_element = matlab.lang.makeValidName(elementName);
    mean_HTT2024_elements.(fieldName_element) = mean(values, 'omitNaN');
end


datasets = {changeFieldNames(mean_HTT2024_elements), changeFieldNames(mean_MG2023_elementalComposition)};  

allFields = {};
for i = 1:length(datasets)
    allFields = [allFields, fieldnames(datasets{i})'];
end

D = regexp(allFields,'[^_]+','once','match');
uniqueFieldsElements = unique(D);

% Calculate the mean values
meanValuesElements = struct(); 

for i = 1:length(uniqueFieldsElements)
    field = uniqueFieldsElements{i};
    values = [];
    for j = 1:length(datasets)
        if isfield(datasets{j}, field)
            values(end+1) = datasets{j}.(field);
        end
    end

    if ~isempty(values)
        meanValuesElements.(field) = mean(values);  % Mean across available datasets
    else
        meanValuesElements.(field) = NaN;  % <-- Optional: fill with NaN if missing everywhere
    end
end


%% Extract CNP As and other elements


originalStructElements  = meanValuesElements;

groupedCompositionStructElements = struct();
groupedCompositionStructElements.Carbon = originalStructElements.Carbon;
groupedCompositionStructElements .Nitrogen = originalStructElements.Nitrogen;
groupedCompositionStructElements .Phosphorous= originalStructElements.Phosphorous;
groupedCompositionStructElements .As  =  originalStructElements.As;
groupedCompositionStructElements .TraceElements  =  0;

structFieldNamesElements   = fieldnames(originalStructElements );
OrganicNames = {'Carbon', 'Nitrogen', 'Phosphorous',  'As'};

for i = 1:length(structFieldNamesElements )
    field = structFieldNamesElements{i};
    if ~ismember(field, OrganicNames)
        value = originalStructElements.(field);
        if ~isnan(value)
            groupedCompositionStructElements.TraceElements = groupedCompositionStructElements.TraceElements + value;
        end  
    end
    
end


ElementsKnown = struct2array(groupedCompositionStructElements); 
ElementsAsh = groupedCompositionStruct.Ash; 
ElementsMoisture =  groupedCompositionStruct.Moisture;
ElementsNotKnown = 100 - sum(sum(ElementsKnown) + ElementsAsh + ElementsMoisture); 


%% Plot Elements
nexttile
dh = bar(categorical({'Elemental Sargassum'}),[ElementsAsh ...
                                                ElementsMoisture ...
                                                ElementsKnown ...
                                                ElementsNotKnown ...
                                                
                                                 ], 'stacked'); 
set(dh, 'FaceColor', 'Flat')

dh(1).CData = min([0.2196    0.0078    0.2314] + 0.4,1) ;  
dh(2).CData = min([0.0081    0.4822    0.6031] + 0.4,1) ; 
for i = 3:length(dh)-1
    dh(i).CData = min(([0.9391    0.3499    0.0103] + (i/10)),1) ;  
end
dh(end).CData = min([0.6377    0.0505    0.1249] + 0.4,1) ;   

ylabel('%biomass gDw')

fields = fieldnames(groupedCompositionStructElements);
fields_clean = strrep(fields, '_', ' ');  % Replace underscores with spaces
newLegendFields = {'Ash','Moisture', 'Carbon', 'Nitrogen', 'Phosphorous',...
    'As', 'Trace elements', 'Organic not known'}; 



lg = legend(newLegendFields, 'Orientation', 'vertical');
    
lg.Layout.Tile = 'north';


%% DEB SARGASSUM
% Obs for any simulation
load('Stylized_facts/GR_Temperatures/obs.mat')
load('Stylized_facts/GR_Temperatures/simu.mat')

obs = obs(3);
simu = simu(3);


finalWd_EC = obs.Wd_EC(end); 
finalWd_EN = obs.Wd_EN(end);
finalWd_V = obs.Wd_V(end); 

final_Wd = [finalWd_EC, finalWd_EN, finalWd_V]; 


perc_gEC = (finalWd_EC / obs.Wd(end))* 100; %  % EC gdW
perc_gEN = (finalWd_EN  / obs.Wd(end)) * 100; %   % EN gdW
perc_gV = (finalWd_V  / obs.Wd(end)) * 100;  %     % MV gdW
%% Plot DEB

nexttile
eh = bar(categorical({'DEB Sargassum'}),[ElementsAsh ...
                                                ElementsMoisture ...
                                                perc_gEC ...
                                                perc_gEN ...
                                                perc_gV ...
                                                 ], 'stacked'); 
set(eh, 'FaceColor', 'Flat')

eh(1).CData = min([0.2196    0.0078    0.2314] + 0.4,1) ;  
eh(2).CData = min([0.0081    0.4822    0.6031] + 0.4,1) ; 
for i = 3:length(eh)
    eh(i).CData = min(([0.9391    0.3499    0.0103] + (i/10)),1) ;  
end

ylabel('%biomass gDw')
newLegendFields = {'Ash','Moisture',"C reserve","N reserve","Structure"};

lg = legend(newLegendFields, 'Orientation', 'vertical');
    
lg.Layout.Tile = 'north';


% legend(newLegendFields); 

%%
set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); set(gcf,'PaperType','A4'); set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')

saveDir = fullfile(pwd,'Figures'); 
figName = 'ElemCompoundComposition.pdf'; 
print(fullfile(saveDir,figName),'-vector','-bestfit','-dpdf')