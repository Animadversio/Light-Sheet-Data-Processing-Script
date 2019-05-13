%% %%%%%%%%%%%%%%%%%%%%%%%%
%%  Delayed Cross Correlation analysis
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

