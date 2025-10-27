Hatt_Lapointe_CNP = [32.4	886 	27.8]; 


Hatt_measured_CNP = [2018	38.54545454545455	2569.53642384106	594.7368421052632;
2019	46.18181818181819	3682.119205298013	742.1052631578948;
2020	62.909090909090914	3443.708609271523	468.42105263157896;
2021	46.18181818181819	2701.9867549668875	515.7894736842105]; 

years = Hatt_measured_CNP(:,1);
CN_measured = Hatt_measured_CNP(:,2);
CP_measured = Hatt_measured_CNP(:,3);
NP_measured = Hatt_measured_CNP(:,4);

% Replicate Hatt_Lapointe_CNP for plotting across the years
CN_ref = repmat(Hatt_Lapointe_CNP(1), size(years));
CP_ref = repmat(Hatt_Lapointe_CNP(2), size(years));
NP_ref = repmat(Hatt_Lapointe_CNP(3), size(years));

% Plot
figure;
subplot(3,1,1);
plot(years, CN_measured, 'o-', 'MarkerSize', 7,...
    'LineWidth', 2, 'Color', [  0.8510    0.3255    0.0980]); 
hold on;
plot(years, CN_ref, '--' , 'LineWidth', 2, 'Color', [0.1294    0.5529    0.6196]);
xticks(years);  % set x-axis ticks to the exact years

ylabel('CN ratio');
ax = gca;  % Find the axes inside the figure
x_limits = xlim(ax);
y_limits = ylim(ax);

xlim(x_limits);
ylim([y_limits(1)- 5  y_limits(2) + 5 ]) ;



subplot(3,1,2);
plot(years, CP_measured, 'o-', 'MarkerSize', 7,...
    'LineWidth', 2, 'Color', [  0.8510    0.3255    0.0980]); 
hold on;
plot(years, CP_ref, '--' , 'LineWidth', 2, 'Color', [0.1294    0.5529    0.6196]);
xticks(years);  % set x-axis ticks to the exact years

ylabel('CP ratio');
ax = gca;  % Find the axes inside the figure
x_limits = xlim(ax);
y_limits = ylim(ax);

xlim(x_limits);
ylim([y_limits(1)- 200  y_limits(2) + 10 ]) ;

subplot(3,1,3);
plot(years, NP_measured, 'o-', 'MarkerSize', 7,...
    'LineWidth', 2, 'Color', [  0.8510    0.3255    0.0980]); 
hold on;
plot(years, NP_ref, '--' , 'LineWidth', 2, 'Color', [0.1294    0.5529    0.6196]);
xticks(years);  % set x-axis ticks to the exact years
ylabel('NP ratio');
ax = gca;  % Find the axes inside the figure
x_limits = xlim(ax);
y_limits = ylim(ax);

xlim(x_limits);
ylim([y_limits(1)- 50  y_limits(2) + 50 ]) ;

legend('Measured ratios Hatt 2024', 'Mean ratios Lapointe 2021', 'Location', 'bestOutside'); 