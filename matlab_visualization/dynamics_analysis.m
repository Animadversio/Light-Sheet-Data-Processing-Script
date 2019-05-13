%% Linear dynamics fitting
% such that  connect_coef * zscore_arr(outperm, t) = 
%           zscore_arr(outperm, t+1) 
connect_coef = ( zscore_arr(outperm, 2:end)  ) / zscore_arr(outperm, 1:end-1); 
%% Calculate Spectrum
e_vec0 = eig(connect_coef(:, :));
%% Spectrum of the Dynamical Matrix 
figure(1);clf;hold on
scatter(real(e_vec0), imag(e_vec0))
ylabel('Imag')
xlabel('Real')
circle(0,0,1);
title('Eigen Spectrum of the 1 order linear dynamical matrix')
axis equal tight
save(1,'Spectrum_1lineards.fig')
saveas(1, 'Spectrum_1lineards.png')
%% visualize dynamics matrix
figure(3)
imagesc(connect_coef,[-0.18,0.18])
colormap('autumn')
colorbar()
axis equal tight
%% coefficient histogram + Spectrum + connection matrix
figure(2);clf
subplot(1,2,1)
imagesc(connect_coef,[-0.18,0.18])
ylabel('Output ROI Index')
xlabel('Input ROI Index')
title('Heatmap of Dynamic Matrix (ROI sorted by linkage)')
colormap('jet')
ch = colorbar();
set(get(ch,'Label'), 'string', 'connection strength (climited)', 'FontSize', 12);
axis equal tight
subplot(1,2,2)
hold on
scatter(real(e_vec0), imag(e_vec0))
ylabel('Imag')
xlabel('Real')
circle(0,0,1);
title('Eigen Spectrum of the 1 order linear dynamical matrix')
axis equal tight
% subplot(1,3,3)
% title('Linear dynamical matrix entry distribution')
% histogram(connect_coef(:))
% xlabel('entry value')
%% Evaluate the fitting (Estimation and Comparison)
%% Assess goodness of fit 
estimated_zscore_arr = connect_coef * zscore_arr(outperm, 1:end-1) ;
estimated_zscore_arr = [zscore_arr(outperm, 1), estimated_zscore_arr];
% err_zscore_arr = estimated_zscore_arr - zscore_arr; 
%% Self autonomous dynamics
autonom_zscore_arr = zeros(size(zscore_arr));
autonom_zscore_arr(:,1) = zscore_arr(outperm,1); 
for t = 2:size(autonom_zscore_arr,2)
    autonom_zscore_arr(:,t) = connect_coef * autonom_zscore_arr(:,t-1);
end
%%
autonom_zscore_arr_o2 = zeros(size(zscore_arr));
autonom_zscore_arr_o2(:,1:2) = zscore_arr(outperm,1:2); 
for t = 3:size(zscore_arr,2)
    autonom_zscore_arr_o2(:, t) = connect_coef_t2 * [autonom_zscore_arr_o2(:,t-2); ...
                                                                                          autonom_zscore_arr_o2(:,t-1)];
end
%%
figure(4)
imagesc(estimated_zscore_arr,[0,5])
%%
figure(5)
imagesc(err_zscore_arr,[0,5])
%%
samp_num = 10;
rand_idx = sort(randperm(size(zscore_arr,1), samp_num)); 
figure(19);clf;hold on
for i=1:samp_num
    subplot(samp_num,1,i);
    hold on;
    plot(zscore_arr(outperm(rand_idx(i)), :), 'LineWidth', 1.5)
    plot(estimated_zscore_arr(rand_idx(i), :), ':', 'LineWidth', 1.5)
    plot(autonom_zscore_arr(rand_idx(i), :), 'k-', 'LineWidth', 1.5)
    plot(autonom_zscore_arr_o2(rand_idx(i), :), 'r-.', 'LineWidth', 1.5)
    ylim([-2,10])
    xlim([0,24000])
    ylabel(num2str(rand_idx(i)))
    if i~=samp_num
        xticklabels([])
    else
        xticks(0:1000:24000)
        xlabel('Time (timestep 0.05s)')
    end
end
suptitle({'Compare Regression and Autonomous Dynamic System Prediction Result (Full Network)',
                'Black trace: Autonomous Prediction; Red -.: 2 order auto prediction', 
                'Orange Dot: Fitting; Blue Trace: original zscore'})





%% Multitimestep regression model
connect_coef_t2 = zscore_arr(outperm, 3:end) / [zscore_arr(outperm, 1:end-2); zscore_arr(outperm, 2:end-1)]; 
%% 
figure(10)
suptitle('Heatmap of Dynamic Matrix of 2 order linear model (ROI sorted by linkage)')
subplot(1,2,1)
imagesc(connect_coef_t2(:, 1:length(outperm)), [-0.18,0.18])%-diag(length(outperm))
ylabel('Output ROI Index')
xlabel('Input ROI Index')
title('Connection Matrix (t-2)')
colormap('jet')
ch = colorbar();
set(get(ch,'Label'), 'string', 'connection strength (climited)', 'FontSize', 12);
axis equal tight
subplot(1,2,2)
imagesc(connect_coef_t2(:, length(outperm)+1:end), [-0.18,0.18])%-diag(length(outperm)
ylabel('Output ROI Index')
xlabel('Input ROI Index')
title('Connection Matrix (t-1)')
colormap('jet')
ch = colorbar();
set(get(ch,'Label'), 'string', 'connection strength (climited)', 'FontSize', 12);
axis equal tight

%% Spectrum of the weight matrix 
e_vec1 = eig(connect_coef_t2(:, 1:length(outperm)));
e_vec2 = eig(connect_coef_t2(:, length(outperm)+1:end));
%%
figure(12)
suptitle('Eigen Spectrum of the 2 order linear dynamical matrix')
subplot(1,2,1)
scatter(real(e_vec1), imag(e_vec1))
ylabel('Imag')
xlabel('Real')
circle(0,0,1);
title('Connection Matrix (t-2)')
axis equal 
subplot(1,2,2)
scatter(real(e_vec2), imag(e_vec2))
ylabel('Imag')
xlabel('Real')
circle(0,0,1);
title('Connection Matrix (t-1)')
axis equal 
%% Spectrum Comparison
figure(13);clf;hold on
suptitle('Eigen Spectrum Comparison of  2 and 1 order linear dynamical matrix')
scatter(real(e_vec1), imag(e_vec1), 16, 'filled', 'MarkerFaceAlpha', 0.75)
scatter(real(e_vec2), imag(e_vec2), 16, 'filled', 'MarkerFaceAlpha', 0.75)
scatter(real(e_vec0), imag(e_vec0), 16, 'filled', 'MarkerFaceAlpha', 0.75)
ylabel('Imag')
xlabel('Real')
circle(0,0,1);
legend({'order2 t-2 matrix', 'order2 t-1 matrix', 'order1 matrix'})
axis equal tight
saveas(13, 'Spectrum_Comparison.png')

%% Local Network Modelling
relind_arr = 717:1313;
loc_connect_coef = ( zscore_arr(outperm(relind_arr), 2:end)  ) / zscore_arr(outperm(relind_arr), 1:end-1); 
loc_connect_coef_t2 = zscore_arr(outperm(relind_arr), 3:end) / [zscore_arr(outperm(relind_arr), 1:end-2); ...
                                                                                                                zscore_arr(outperm(relind_arr), 2:end-1)]; 
loc_e_vec0 = eig(loc_connect_coef(:, :));
loc_e_vec1 = eig(loc_connect_coef_t2(:, 1:length(relind_arr)));
loc_e_vec2 = eig(loc_connect_coef_t2(:, length(relind_arr)+1:end));
%%
figure(17);clf;hold on
title({'Eigen Spectrum Comparison of  2 and 1 order linear dynamical matrix', ...
    '(outperm(717:1313) subnetwork)'})
scatter(real(loc_e_vec1), imag(loc_e_vec1), 16, 'filled', 'MarkerFaceAlpha', 0.9)
scatter(real(loc_e_vec2), imag(loc_e_vec2), 16, 'filled', 'MarkerFaceAlpha', 0.9)
scatter(real(loc_e_vec0), imag(loc_e_vec0), 16, 'filled', 'MarkerFaceAlpha', 0.9)
ylabel('Imag')
xlabel('Real')
circle(0,0,1);
legend({'order2 t-2 matrix', 'order2 t-1 matrix', 'order1 matrix'})
axis equal tight
saveas(gcf, 'Loc717-1313_Spectrum_Comparison.png')
%%
figure()
suptitle({'Heatmap of Dynamic Matrix of 2 order linear model (ROI sorted by linkage)',...
                '(outperm(717:1313) subnetwork)'})
subplot(1,3,1)
imagesc(relind_arr, relind_arr, loc_connect_coef(:, 1:length(relind_arr)), [-0.18,0.18])%-diag(length(outperm))
ylabel('Output ROI Index')
xlabel('Input ROI Index')
title('Connection Matrix (order 1)')
colormap('jet')
ch = colorbar();
set(get(ch,'Label'), 'string', 'connection strength (climited)', 'FontSize', 12);
axis equal tight
subplot(1,3,2)
imagesc(relind_arr, relind_arr, loc_connect_coef_t2(:, 1:length(relind_arr)), [-0.18,0.18])%-diag(length(outperm))
ylabel('Output ROI Index')
xlabel('Input ROI Index')
title('Connection Matrix (t-2)')
colormap('jet')
ch = colorbar();
set(get(ch,'Label'), 'string', 'connection strength (climited)', 'FontSize', 12);
axis equal tight
subplot(1,3,3)
imagesc(relind_arr, relind_arr, loc_connect_coef_t2(:, length(relind_arr)+1:end), [-0.18,0.18])%-diag(length(outperm)
ylabel('Output ROI Index')
xlabel('Input ROI Index')
title('Connection Matrix (t-1)')
colormap('jet')
ch = colorbar();
set(get(ch,'Label'), 'string', 'connection strength (climited)', 'FontSize', 12);
axis equal tight
%% Linear 1 step estimation 
loc_est_act_arr = zeros(size(loc_auto_act_arr));
loc_est_act_arr(:, 1) = zscore_arr(outperm(relind_arr), 1);
loc_est_act_arr(:,2:end) = loc_connect_coef * zscore_arr(outperm(relind_arr), 1:end-1);
%% Autonomous Dynamic system
loc_auto_act_arr = zeros(length(relind_arr), size(zscore_arr,2));
loc_auto_act_arr(:, 1) = zscore_arr(outperm(relind_arr), 1);
for  t = 2:size(zscore_arr,2)
    loc_auto_act_arr(:,t) = loc_connect_coef * loc_auto_act_arr(:,t-1); 
end
%% 2 order prediction
loc_auto_act_arr_o2 = zeros(length(relind_arr), size(zscore_arr,2));
loc_auto_act_arr_02(:, 1:2) = zscore_arr(outperm(relind_arr), 1:2);
for  t = 3:size(zscore_arr,2)
    loc_auto_act_arr_02(:,t) = loc_connect_coef_t2* [loc_auto_act_arr_02(:,t-2); ...
                                                                                            loc_auto_act_arr_02(:,t-1)];
end
%%
samp_num = 10;
rand_idx = sort(randperm(length(relind_arr), samp_num)); %randi(length(rand_idx),[1,10]); 
figure(5);clf;hold on
for i=1:samp_num
    subplot(samp_num,1,i);
    hold on;
    plot(zscore_arr(outperm(relind_arr(rand_idx(i))), :), 'LineWidth', 1.5)
    plot(loc_est_act_arr(rand_idx(i), :), ':', 'LineWidth', 1.5)
    plot(loc_auto_act_arr(rand_idx(i), :),'k-', 'LineWidth', 1.5)
    plot(loc_auto_act_arr_02(rand_idx(i), :),'r-.', 'LineWidth', 1.5)
    ylim([-2,10])
    xlim([0,24000])
    ylabel(num2str(relind_arr(rand_idx(i))))
    if i~=samp_num
        xticklabels([])
    else
        xticks(0:1000:24000)
        xlabel('Time (timestep 0.05s)')
    end
end
suptitle({'Compare Regression and Autonomous Dynamic System Prediction Result (717-1313 subnet)',
                'Black trace: 1 order Autonomous Prediction; Red -.: 2 order auto prediction',
                'Orange Dot: Fitting; Blue Trace: original zscore'})

%% Auto regressive model 

%% L1 penalizaed regression model 
connect_coef_l1 = lassoglm(zscore_arr(outperm, 1:end-1),zscore_arr(outperm, 2:end));

