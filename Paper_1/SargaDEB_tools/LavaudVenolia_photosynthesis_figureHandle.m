data_I_J_O2_Venolia = [ ...  % Light intensity [mol quanta m-2 h-1], gross photosynthesis rate [mg O2 gdW-2 h-1]  %Data obtained from Roman
0.03225806451612903 0.010067114093959731
0.025806451612903226 0.8053691275167785
0.06451612903225806 1.308724832214765
0.06451612903225806 2.0738255033557045
0.17419354838709677 2.6275167785234896
0.17419354838709677 3.483221476510067
0.3548387096774194 3.6845637583892614
0.3548387096774194 4.419463087248322
1.0774193548387097 4.0369127516778525
1.0774193548387097 4.600671140939597
2.161290322580645 4.4597315436241605
2.161290322580645 5.1342281879194624
3.238709677419355 4.640939597315436
3.232258064516129 5.023489932885906
3.232258064516129 5.476510067114093];

data_I_J_O2_Venolia_model = [ ...
    0.03225806451612903	 0.010067114093959731
0.025806451612903226	 0.8053691275167785
0.06451612903225806	 1.308724832214765
0.13548387096774195	 2.0838926174496644
0.2064516129032258	 2.657718120805369
0.3935483870967742	 3.4530201342281877
0.7290322580645161	 4.12751677852349
1.335483870967742	 4.580536912751677
1.9161290322580644	 4.771812080536913
2.7096774193548385	 4.912751677852349 
3.238709677419355	 4.9630872483221475]; 





Lavaud_25_data = [
1.8495684340320593, 0.009404388714733543
49.9383477188656, 0.7711598746081505
86.92971639950679, 0.9404388714733543
109.1245376078915, 1.3636363636363638
131.3193588162762, 1.8056426332288402
175.70900123304563, 1.9184952978056427
262.6387176325524, 2.0595611285266457
442.04685573366214, 2.1912225705329154
619.6054254007398, 1.8996865203761757
811.960542540074, 1.7210031347962382
1074.5992601726264, 1.64576802507837
1222.5647348951911, 1.4670846394984327
1414.9198520345253, 1.090909090909091];


Lavaud_25_model = [
0, 0.018808777429467086
7.398273736128237, 0.4137931034482759
11.097410604192355, 0.8369905956112853
20.34525277435265, 1.3260188087774296
73.98273736128237, 1.6175548589341693
159.0628853267571, 1.7115987460815048
266.3378545006165, 1.7210031347962382
369.9136868064119, 1.730407523510972
467.940813810111, 1.749216300940439
604.8088779284834, 1.7398119122257054
712.0838471023428, 1.730407523510972
830.4562268803946, 1.730407523510972
950.6781750924785, 1.749216300940439
1052.4044389642418, 1.7398119122257054
1167.0776818742295, 1.7398119122257054
1298.3970406905057, 1.749216300940439
1355.7336621454995, 1.749216300940439
1398.2737361282368, 1.7586206896551726];


% Conversion: multiply all y-values (column 2) by 31.25
data_I_J_O2_Venolia(:,2) = data_I_J_O2_Venolia(:,2) * 31.25; % Convert mg02 into micromol 02
data_I_J_O2_Venolia_model(:,2) = data_I_J_O2_Venolia_model(:,2) * 31.25;  %Convert mg02 into micromol 02

Lavaud_25_data(:,2) = Lavaud_25_data(:,2) * 31.25;  %Convert mg02 into micromol 02
Lavaud_25_model(:,2) = Lavaud_25_model(:,2) * 31.25;  %Convert mg02 into micromol 02

% Split the data into x and y
x = data_I_J_O2_Venolia_model(:,1) / 1e-6 / 3600; %convert to micro and to per second
y = data_I_J_O2_Venolia_model(:,2);

xL = Lavaud_25_model(:,1) ; 
yL = Lavaud_25_model(:,2) ; 

% Create a finer x-axis for interpolation
xq = linspace(min(x), max(x), 300);  % 300 points for smoothness
xqL = linspace(min(xL), max(xL), 300);  % 300 points for smoothness

% Interpolate using spline (smooth) or 'pchip' (shape-preserving)
yq = interp1(x, y, xq, 'pchip');
yq(:,1) = 0.010067114093959731; 

yqL = interp1(xL, yL, xqL, 'pchip');
yqL(:,1) = 0.010067114093959731; 


% Plot
figure(10);
plot(xq, yq, '-', 'LineWidth', 2, 'Color', [0.0745    0.6235    1.0000]);      hold on;  % Interpolated smooth curve
scatter(data_I_J_O2_Venolia(:,1) / 1e-6 / 3600, data_I_J_O2_Venolia(:,2),  'MarkerEdgecolor', 'none', ...
    'MarkerFaceColor',  [0.0745    0.6235    1.0000]);


plot(xqL, yqL, '-', 'LineWidth', 2, 'Color', [1.0000    0.4118    0.1608]);      hold on;  % Interpolated smooth curve
scatter(Lavaud_25_data (:,1), Lavaud_25_data (:,2),  'MarkerEdgecolor', 'none', ...
    'MarkerFaceColor',  [1.0000    0.4118    0.1608]);

% legend('Parameter estimation ', 'S. fluitans at 23°C HL acclimated','Sacharinna at 13°C','', 'Ulva at 25°C', '')
xlabel('micro mol quanta m-2 s-1')
ylabel('micro mol O2 gdW-1 h-1')
% xlim(x_limits);
% ylim(y_limits);
%%
% r = figure(5);
% ax = findall(r, 'Type', 'axes');  % Find the axes inside the figure
% x_limits = xlim(ax);
% y_limits = ylim(ax);
% 
% xlim(x_limits);
% ylim(y_limits);
%%
f = figure(4);                % use the correct figure number
ax = findall(f, 'Type', 'axes');
xlim(ax, x_limits);
ylim(ax, y_limits);
% find all plotted objects in the axes
lines = findall(ax, 'Type', 'line');  % or 'Type', 'patch', etc. if relevant

% set legend (the order might be reversed, so double-check)
% legend(ax,  [lines(2), lines(4)], {'Low light acclimated', 'High light acclimated', ''}); 

%%
%Estimate Irradiance as a function of light intensity and t based on
%Fourier equations

t = 1:24; 

lightIntensity = 60 * 10^6 / 86400; %60 ein m-2 d-1, 10e6 to convert into micro mol and 86400 seconds in one day


Irradiance = get_irradiance(t, lightIntensity) ; 
irradianceday = Irradiance(:,9:18); 


LP1986 = [... %hour of the day , mg C g Wd-1 h-1
    896.830985915493, 0.6818980667838312
996.1267605633802, 0.8646748681898067
1096.4788732394366, 0.9138840070298769
1189.4366197183099, 0.9490333919156414
1296.1267605633802, 1.0298769771528997
1391.1971830985915, 0.9103690685413005
1497.887323943662, 0.7627416520210896
1597.1830985915492, 0.59402460456942
1699.6478873239435, 0.17574692442882248
1802.112676056338, 0.0984182776801406]; 


LP1986_micromol = LP1986(:,2)' .* 31.25; % Convert mg02 into micromol 02

plot(1:24, Irradiance)


plot(Irradiance(9:18), LP1986_micromol)

scatter(LP1986(:,1), LP1986_micromol)