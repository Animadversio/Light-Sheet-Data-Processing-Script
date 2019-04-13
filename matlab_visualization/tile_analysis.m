%% Histogram / scatter of statistical features across tiles
histogram2(stats_arr(:,1),stats_arr(:,3),50)
%%
figure(5)
histogram(log(stats_arr(:,11)),50)
%%
raw_traces = raw_traces';
%%
figure()
plot(raw_traces(967,:))
%% Show the activity traces of all the tiles 
imagesc(zscore_arr,[0,5]) % crop the clim to 5 std
title("zscore activity trace")
ylabel("Tile id")
xlabel("Time")
%%
figure(3)
plt = plot(zscore_arr(idx,7035:7445)');
%alpha(.5)

%% Moving average the zscore trace 
mean_arr = movmean(zscore_arr, 11, 2, 'Endpoints', 'shrink');
pop_activity = sum(mean_arr, 1); % summed population activity 
pop_active_unit = sum(mean_arr > 5, 1); % unit number of significant activity 
%% 
[pks,locs] = findpeaks(pop_active_unit, 'MinPeakDistance', 30, ...
        'MinPeakProminence', 10,'Annotate','extents');
%% Findpeaks in population data and show the traces! 
figure(7);clf;hold on
subplot(2,1,1);plot(pop_activity);ylabel("Summed population acitvity");
subplot(2,1,2);hold on
plot(pop_active_unit);ylabel("Thresholded active unit #");
scatter(locs, pks,'^')
text(locs, pks, arrayfun(@(l) sprintf('%d', l),  ...
        pks, 'UniformOutput', false))
% output = cellfun(@(l) sprintf('Peak = %d', l), pks, 'UniformOutput', false

%% For each peak in population trace, extract the activity clip around it 
%  for all the active units 
for i=1:numel(locs)
    pk_loc = uint16(locs(i));
    activ_num = pks(i);
    tile_idx = uint16(find(mean(mean_arr(:,pk_loc-5:pk_loc+5),2) > 5));
    window = max(1, pk_loc-200):min(pk_loc+500,size(zscore_arr, 2));
    h = figure(8);set(h,'position',[34         107        1407         691]);
    plot(window, zscore_arr(tile_idx, window)')
    title(sprintf("Active trace around %d (Unit # %d)", pk_loc, activ_num))
    xlabel("Time point");ylabel("z score")
    xlim([window(1), window(end)+1])
    YL = ylim();
    ylim([0, YL(2)])
    saveas(h, sprintf('/Users/binxu/HolyLab/Synfire_peak_%05d.png',pk_loc))
end
%% Lowpass filter trace
trace = zscore_arr(tile_idx(20), :);
y = filter1('lp',double(trace), 'fs', 20, 'fc', 0.5);
figure(4);clf;hold on;
plot(trace);
plot(y)
%% Save the lowpass filter traces 
lpf_zscore_arr = zeros(size(zscore_arr));
for i=1:size(zscore_arr)
    trace = zscore_arr(i,:);
    lpf_zscore_arr(i,:) = filter1('lp',double(trace), 'fs', 20, 'fc', 1);
end
for i=1:numel(locs)
    pk_loc = uint16(locs(i));
    activ_num = pks(i);
    tile_idx = uint16(find(mean(lpf_zscore_arr(:,pk_loc-5:pk_loc+5),2) > 5));
    window = max(1, pk_loc-200):min(pk_loc+500,size(zscore_arr, 2));
    h = figure(8);set(h,'position',[34         107        1407         691]);
    plot(window, lpf_zscore_arr(tile_idx, window)')
    title(sprintf("Active trace around %d (Unit # %d)", pk_loc, activ_num))
    xlabel("Time point");ylabel("z score")
    xlim([window(1), window(end)+1])
    YL = ylim();
    ylim([0, YL(2)])
    saveas(h, sprintf('/Users/binxu/HolyLab/Synfire_peak_lpf_%05d.png',pk_loc))
end
clear lpf_zscore_arr
%%
% y = lowpass(zscore_arr(781,:),5,20);
% figure(9);hold on
% plot(zscore_arr(781,:))
% plot(y)
%% Draw Correlogram
corr_mat = corrcoef(zscore_arr');
%% 
figure()
imagesc(corr_mat)
%% Hirarchical Clustering and 
% Y = pdist(zscore_arr,'correlation');
% Z = linkage(Y,'average');
lpf_zscore_arr = zeros(size(zscore_arr)); % Lowpass filtered version of zscore trace
for i=1:size(zscore_arr)
    trace = zscore_arr(i,:);
    lpf_zscore_arr(i,:) = filter1('lp',double(trace), 'fs', 20, 'fc', 1);
end
Z = linkage(lpf_zscore_arr,'average','correlation');
% Z = linkage(zscore_arr,'average','correlation');  % Another version is to
% use the unfiltered trace ! 
%% Inspect the dendrogram of the linkaged data 
figure(12)
dendrogram(Z,150)
%%
T = cluster(Z,'Cutoff',0.5,'Criterion','distance');
% T = clusterdata(zscore_arr,'Criterion','distance','Cutoff',0.5);
[sort_class, perm_ind] = sort(T);
fprintf("Cluster number %d\n",max(T))
%% Count the tile numbers belonging to a cluster
C = unique(sort_class);
C_cnt=zeros(1,length(C));
for i = 1:length(C)
    C_cnt(i) = sum(sort_class==C(i));
end
[~, perm_Cid] = sort(C_cnt, 'descend');
C = C(perm_Cid);
C_cnt = C_cnt(perm_Cid); % re arrange the Classes from the largest cluster to the smallest 
fprintf("Maximum cluster sizes %s\n",num2str(C_cnt(1:5)))
%% Show the re-arranged traces 
figure(11)
imagesc(corr_mat(perm_ind, perm_ind))

%% 
C_id = 1;
cluster_tile_id = int16(find(T==C(C_id)));
figure(10);
plot(lpf_zscore_arr(cluster_tile_id,:)')
title(sprintf('Cluster %d (tile number %d)',C(C_id),C_cnt(C_id)))
%% spatial location of the tiles 

% C_id = 5;
C_id_list = [2,3,4,5,6,7,8,9,10,11,12,13,14];
C_color_list = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
    [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], ...
    [0.6350 0.0780 0.1840], 'red', 'green', 'blue', 'magenta', 'cyan', 'yellow'};
for i = 1:numel(C_id_list)
C_id = C_id_list(i);
cluster_tile_id = int16(find(T==C(C_id)));
% cluster_tile_id = int16([967, 968, 969, 1000]);
vol_canvas = zeros(400,198,40);
for tid = cluster_tile_id'
    global_tid = tile_idx{tid};% tid;%
%     disp(global_tid)
    pos_rng = tile_pos{global_tid};
    vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
            pos_rng{2}(1):pos_rng{2}(2), ...
            pos_rng{3}(1):pos_rng{3}(2)) = tile_space{global_tid} + ...
            vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
            pos_rng{2}(1):pos_rng{2}(2), ...
            pos_rng{3}(1):pos_rng{3}(2));
end
%
% pcolor3(vol_canvas)
%
figure(1);histogram(vol_canvas(vol_canvas>0))
%
rX = (1:size(vol_canvas,1))*0.65;
rY = (1:size(vol_canvas,2))*0.65;
rZ = (1:size(vol_canvas,3))*5;
figure(2);hold on 
p = patch(isosurface(rY,rX,rZ,vol_canvas,0.035));
isonormals(rY,rX,rZ,vol_canvas,p)
p.FaceColor = C_color_list{i};%'magenta';
p.EdgeColor = 'none';
p.FaceAlpha = 0.4;
end
%
xlabel("Y")
ylabel("X")
zlabel("Z")
daspect([1 1 1])
view(3); 
axis tight

