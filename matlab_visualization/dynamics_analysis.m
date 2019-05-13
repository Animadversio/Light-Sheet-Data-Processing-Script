%% Linear dynamics fitting
connect_coef = ( zscore_arr(outperm, 2:end)  ) / zscore_arr(outperm, 1:end-1); 
% such that  connect_coef * zscore_arr(outperm, t) = 
%           zscore_arr(outperm, t+1) - zscore_arr(outperm, t)
%% Assess goodness of fit 
estimated_zscore_arr = connect_coef * zscore_arr(outperm, 1:end-1) ;
estimated_zscore_arr = [zscore_arr(outperm, 1), estimated_zscore_arr];
err_zscore_arr = estimated_zscore_arr - zscore_arr; 
%% Spaectrum of the Dynamical Matrix 
figure()
e_vec0 = eig(connect_coef(:, :));
scatter(real(e_vec0), imag(e_vec0))
ylabel('Imag')
xlabel('Real')
title('Spectrum of the 1 order linear dynamical matrix')
axis equal 
%%
figure(4)
imagesc(estimated_zscore_arr,[0,5])
%%
figure(5)
imagesc(err_zscore_arr,[0,5])
%%
rand_idx = randi(2271,[1,10]); 
figure(6);clf;hold on
for i=1:numel(rand_idx)
    subplot(numel(rand_idx),1,i);
    hold on;
    plot(estimated_zscore_arr(rand_idx(i), :))
    plot(zscore_arr(outperm(rand_idx(i)), :))
    plot(autonom_zscore_arr(rand_idx(i), :))
    ylim([-2,10])
end
%% Self autonomous dynamics
autonom_zscore_arr = zeros(size(zscore_arr));
autonom_zscore_arr(:,1) = zscore_arr(outperm,1); 
for t = 2:size(autonom_zscore_arr,2)
    autonom_zscore_arr(:,t) = connect_coef * autonom_zscore_arr(:,t-1);
end

%% visualize dynamics 
figure(3)
imagesc(connect_coef,[-0.18,0.18])
colormap('autumn')
colorbar()
axis equal tight
%% coefficient histogram
figure()
histogram(connect_coef(:))

%% L1 penalizaed regression model 
connect_coef_l1 = lassoglm(zscore_arr(outperm, 1:end-1),zscore_arr(outperm, 2:end));
%% Multitimestep regression

connect_coef_t2 = zscore_arr(outperm, 3:end) / [zscore_arr(outperm, 1:end-2); zscore_arr(outperm, 2:end-1)]; 


%% 
figure(10)
subplot(1,2,1)
imagesc(connect_coef_t2(:, 1:length(outperm)), [-0.01,0.01])%-diag(length(outperm))
axis equal tight
colorbar()
subplot(1,2,2)
imagesc(connect_coef_t2(:, length(outperm)+1:end), [-0.01,0.01])%-diag(length(outperm)
axis equal tight
colorbar()
figure(11)
subplot(1,2,1)
imagesc(connect_coef_t2(:, 1:length(outperm)), [-0.01,0.01])%-diag(length(outperm))
axis equal tight
colorbar()
subplot(1,2,2)
imagesc(connect_coef_t2(:, 1:length(outperm))', [-0.01,0.01])%-diag(length(outperm)
axis equal tight
colorbar()
%% Spectrum of the weight matrix 

figure(12)
subplot(1,2,1)
e_vec = eig(connect_coef_t2(:, 1:length(outperm)));
scatter(real(e_vec), imag(e_vec))
axis equal 
subplot(1,2,2)
e_vec2 = eig(connect_coef_t2(:, length(outperm)+1:end));
scatter(real(e_vec2), imag(e_vec2))
axis equal 
%% Auto regressive model 


%%  delayed correlation 
tic
lag_mat = zeros(size(zscore_arr,1),'int16'); % uint8
corr_ef_mat = zeros(size(zscore_arr,1),'single');
for i = 1:size(zscore_arr,1)
    disp(i)
    parfor j = i+1:size(zscore_arr,1)
        [corr,lag] = xcorr(zscore_arr(i, :)',zscore_arr(j, :)',30,'coeff' );
        [m_cf, mid] = max(corr);
        corr_ef_mat(i,j) = m_cf;
        lag_mat(i,j) = lag(mid);
%         corr_ef_mat(j,i) = m_cf;
%         lag_mat(j,i) = - lag(mid);
    end
end 
toc
figure(8)
imagesc(lag_mat)
figure(9)
imagesc(corr_ef_mat)
save('xcorr_mat.mat', 'lag_mat', 'corr_ef_mat')
%%
lag_mat_sym = int16(lag_mat) - int16(lag_mat');
corr_ef_mat_sym = corr_ef_mat + corr_ef_mat';
%%
figure(8)
imagesc(lag_mat_sym(outperm,outperm))
figure(9)
imagesc(corr_ef_mat_sym(outperm,outperm))
%%
tic
idx_array = outperm(717:1313);%(2007:2040);
lag_submat = zeros(length(idx_array),'int16'); % uint8
corr_ef_submat = zeros(length(idx_array),'single');
for i = 1:length(idx_array)
    disp(i)
    parfor j = i+1:length(idx_array)
        [corr,lag] = xcorr(zscore_arr(idx_array(i), :)',zscore_arr(idx_array(j), :)',50,'coeff' );
        [m_cf, mid] = max(corr);
        corr_ef_submat(i,j) = m_cf;
        lag_submat(i,j) = lag(mid);
%         corr_ef_submat(j,i) = m_cf;
%         lag_submat(j,i) = - lag(mid);
    end
end 
toc
corr_ef_submat = corr_ef_submat + corr_ef_submat';
lag_submat = lag_submat - lag_submat';
figure(11);clf;
imagesc(corr_ef_submat)
colorbar()
axis equal tight
figure(12);clf;
imagesc(lag_submat)
colorbar()
axis equal tight
save('xcorr_submatrix_717-1313.mat','corr_ef_submat', 'lag_submat')
% save('xcorr_submatrix_2007-2040.mat','corr_ef_submat', 'lag_submat')

figure(16);clf
subplot(2,2,1)
plot(zscore_arr(idx_array,:)')
xlabel('time')
ylabel('zscore F')
xlim([0,24000])
ylim([-0.1,20])
subplot(2,2,3)
imagesc(zscore_arr(idx_array,:))
xlabel('time')
ylabel('ROI id')
subplot(2,2,2)
imagesc(corr_ef_submat)
ch1 = colorbar();
set(get(ch1,'Label'), 'string', 'maximum cross correlation', 'FontSize', 12);
axis equal tight
subplot(2,2,4)
imagesc(lag_submat)
ch2 = colorbar();
set(get(ch2,'Label'), 'string',  'lag time for maximum correlation', 'FontSize', 12);
axis equal tight
suptitle("Cluster Cross Correlation Plot (outperm(717-1313))")

