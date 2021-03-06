%% Import data
if ismac
    cd '/Users/binxu/Library/Mobile Documents/com~apple~CloudDocs/Data_Transport/'
elseif ispc
    cd 'C:\Users\binxu\iCloudDrive\Data_Transport'
end
load('raw_traces.mat')
load('tile_pos2.mat')
load('tile_space.mat')
load('tile_traces.mat')
if ismac
    cd '/Users/binxu/Holy_Optical_Imaging/matlab_visualization'
elseif ispc
    cd 'D:\Light-Sheet-Data-Processing-Script\matlab_visualization'
end
load('image_vol.mat')
tile_idx = cell2mat(tile_idx); 
raw_traces = raw_traces';
%%
figure()
imagesc(raw_traces, [0,0.02])
xticks(0:1000:24000)
ylabel('Tile id ', 'FontSize', 14)
xlabel('Time (timestep 0.05s)', 'FontSize', 14)
title('Raw activity trace ', 'FontSize', 20)
ch = colorbar();
set(get(ch,'Label'), 'string', 'raw activity limited', 'FontSize', 14);
%% Histogram / scatter of statistical features across tiles
histogram2(stats_arr(:,1),stats_arr(:,3),50)
%%
figure(5)
histogram(log(stats_arr(:,11)),50)
%%
figure(7)
imagesc(raw_traces(tile_idx(outperm), :), [0,0.02])
xticks(0:1000:24000)
ylabel('Tile id (sorted by linkage)', 'FontSize', 14)
xlabel('Time (timestep 0.05s)', 'FontSize', 14)
title('raw activity trace ', 'FontSize', 20)
ch = colorbar();
set(get(ch,'Label'), 'string', 'raw activity limited', 'FontSize', 14);
%%
figure()
plot(raw_traces(967,:))
%% Show the activity traces of all the tiles 
figure(2)
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
%% Calculate Correlogram
corr_mat = corrcoef(zscore_arr');
%% 
figure(1)
imagesc(corr_mat)
title("Correlogram of all selected tiles",'fontsize',18);
set(gcf, 'position', [100, 0, 1022, 797])
colorbar()
axis equal tight
saveas("correlogram.png")
%% Hirarchical Clustering and sorting of (lpf_zscore_arr)
lpf_zscore_arr = zeros(size(zscore_arr)); % Lowpass filtered version of zscore trace
for i=1:size(zscore_arr)
    trace = zscore_arr(i,:);
    lpf_zscore_arr(i,:) = filter1('lp',double(trace), 'fs', 20, 'fc', 1);
end
Z = linkage(lpf_zscore_arr,'average','correlation');
% Y = pdist(zscore_arr,'correlation');
% Z = linkage(Y,'average');
% Z = linkage(zscore_arr,'average','correlation');  % Another version is to
% use the unfiltered trace ! 
%% Inspect the dendrogram of the linkaged data 
figure(15)
[~, T, outperm] = dendrogram(Z,0); 
title("Dendrogram of (low-pass filtered traces, correlation, average link) linkage", 'FontSize', 18)
ylabel("Link height", 'FontSize', 18)
set(gcf, 'position', [34         360        1407         438])
%% Visualize
figure(24)
set(gcf, 'position', [1          41        1920         963])
subplot('Position', [0.05, 0.7, 0.9, 0.2])
[~, T, outperm] = dendrogram(Z,0, 'ColorThreshold', 0.5); % , 'Orientation', 'left'
title(["Dendrogram", "(low-pass filtered z-score traces, correlation, average linkage)"], 'FontSize', 18)
ylabel("Link height", 'FontSize', 18)
xticklabels([])
xticks([])
subplot('Position', [0.05, 0.08, 0.9, 0.6])
imagesc(zscore_arr(outperm,:)',[0,5])
xlabel('Sorted Tile id', 'FontSize', 18)
ylabel('Time', 'FontSize', 18)
% saveas(24, 'Dendrogram_sorted_trace.png')
%%
corr_mat = corrcoef(zscore_arr');
%%
% get the index sequence from the full dendrogram tree and used to sort the
figure()
imagesc(corr_mat(outperm, outperm))
title("Correogram of all selected tiles (sorted by linkage)",'fontsize',18);
set(gcf, 'position', [100, 0, 1022, 797])
ch = colorbar();
set(get(ch,'Label'), 'string', 'correlation coefficient', 'FontSize', 12);
axis equal tight
%%
figure(6)
imagesc(zscore_arr(outperm,:),[0,5]) % crop the clim to 5 std
title("zscore activity trace (sorted by linkage)",'FontSize',18)
ylabel("Tile id")
xlabel("Time")
set(17,'position',[34,0,1400,800])
%%
% hline(separator_h,'w:')
%% Cutofff the tree to get cluster
T = cluster(Z,'Cutoff',0.5,'Criterion','distance');
% T = clusterdata(zscore_arr,'Criterion','distance','Cutoff',0.5);
[sort_class, perm_ind] = sort(T);
fprintf("Cluster number %d\n",max(T))
%% Show the re-arranged traces 
figure(11)
imagesc(corr_mat(perm_ind, perm_ind))
title("Correogram of all selected tiles (clustered)",'fontsize',18);
set(gcf, 'position', [100, 0, 1022, 797])
colorbar()
axis equal tight
saveas(11,"correlogram_rearrange.png")
%% Count the tile numbers belonging to a cluster
C = unique(sort_class);
C_cnt = zeros(1,length(C));
for i = 1:length(C)
    C_cnt(i) = sum(sort_class==C(i));
end
[~, perm_Cid] = sort(C_cnt, 'descend');
C = C(perm_Cid);
C_cnt = C_cnt(perm_Cid); % re arrange the Classes from the largest cluster to the smallest 
fprintf("Maximum cluster sizes %s\n",num2str(C_cnt(1:5)))

%% Sort the clusters by the size of cluster
perm_ind_sort = [];
separator_h = [];
for cid = C'
    perm_ind_sort = [perm_ind_sort;perm_ind(sort_class==cid)];
    separator_h = [separator_h,length(perm_ind_sort)+0.5];
end
%% Show the heatmap re-arranged by clusters
figure(3)
imagesc(zscore_arr(perm_ind_sort,:),[0,5]) % crop the clim to 5 std
title("zscore activity trace (sorted by cluster size)",'FontSize',18)
ylabel("Tile id")
xlabel("Time")
hline(separator_h,'w:')
%% 
figure(5)
imagesc(corr_mat(perm_ind_sort, perm_ind_sort))
title("Correlogram of all selected tiles (sorted by cluster size)",'fontsize',18);
set(gcf, 'position', [100, 0, 1022, 797])
colorbar()
axis equal tight
saveas(5,"correlogram_sort_size.png")
saveas(5,"correlogram_sort_size.fig")
%% plot all the traces belonging to a cluster
C_id = 4;
cluster_tile_id = int16(find(T==C(C_id)));
figure(10);
plot(lpf_zscore_arr(cluster_tile_id,:)')
title(sprintf('Cluster %d (tile number %d)',C(C_id),C_cnt(C_id)))
%% Interactively inspect each peak! 
ind_array = 1025:1196;
figure(1);
plot(lpf_zscore_arr(outperm(1025:1196),:)','-','LineWidth',0.5)
xlim([0, 24000])
xticks()
xlabel("Time (timestep 0.05s)")
ylabel("zscore (low pass filtered)")
title(sprintf('Cluster %d-%d', ind_array(1), ind_array(end)))
%% Isosurface Visualization of Spatial Location of the tiles 
% C_id = 5;
C_id_list = int16([2,]);%[2,3,4,5,6,7,8,9,10,11,12,13,14];
C_color_list = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
    [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], ...
    [0.6350 0.0780 0.1840], 'red', 'green', 'blue', 'magenta', 'cyan', 'yellow'};
for i = 1:numel(C_id_list)
C_id = C_id_list(i);
cluster_tile_id = int16(find(T==C(C_id)));
% cluster_tile_id = int16([967, 968, 969, 1000]);
vol_canvas = zeros(400,198,40);
for tid = cluster_tile_id'
    global_tid = tile_idx{tid};% tid;
%     disp(global_tid)
    pos_rng = tile_pos{global_tid};
    vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
            pos_rng{2}(1):pos_rng{2}(2), ...
            pos_rng{3}(1):pos_rng{3}(2)) = tile_space{global_tid} + ...
            vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
            pos_rng{2}(1):pos_rng{2}(2), ...
            pos_rng{3}(1):pos_rng{3}(2));
end
% pcolor3(vol_canvas)
figure(2);histogram(vol_canvas(vol_canvas>0)) % show the histogram of volume values
%
rX = (1:size(vol_canvas,1))*0.65;
rY = (1:size(vol_canvas,2))*0.65;
rZ = (1:size(vol_canvas,3))*5;
figure(3);hold on 
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
%% Overlay image volume and ROI tiles 
%% Load the image file and import the anatomical volumes
m = memmapfile('/Volumes/Seagate_Backup_Binxu/tmp_tfm_img3.bin','Format',{'single',[400, 198, 40, 24000],'x'});
figure(5);imagesc(m.Data.x(:,:,10,1),'gray');axis equal tight
image_vol = m.Data.x(:,:,:,1800);
%% Load the saved image volume
load('image_vol.mat')
%% Pre-processing and limit the range of values
UL = 0.008;BL = 0.0015;
image_vol_crop = image_vol;
image_vol_crop(image_vol_crop>UL) = UL;
image_vol_crop(image_vol_crop<BL) = 0;
%% Make the label canvas and the volume 
thresh = 0.05;
label_canvas = zeros(400,198,40,'int16');
vol_canvas = zeros(400,198,40,'double');
C_id_list = int16(2:30);
for i = 1:numel(C_id_list)
    C_id = C_id_list(i);
    cluster_tile_id = int16(find(T==C(C_id)));
    for tid = cluster_tile_id'
        global_tid = tile_idx{tid}; 
        tile_pad = C_id*int16(tile_space{global_tid} > thresh);
        pos_rng = tile_pos{global_tid}; 
        label_canvas(pos_rng{1}(1):pos_rng{1}(2), ... 
                pos_rng{2}(1):pos_rng{2}(2), ... 
                pos_rng{3}(1):pos_rng{3}(2)) = tile_pad;%C_id;
        vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
                pos_rng{2}(1):pos_rng{2}(2), ...
                pos_rng{3}(1):pos_rng{3}(2)) = tile_space{global_tid} + ...
                vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
                pos_rng{2}(1):pos_rng{2}(2), ...
                pos_rng{3}(1):pos_rng{3}(2));
    end
end
% volshow(vol_canvas,label_canvas,'ScaleFactors',[0.65,0.65,5])
volumeViewer(image_vol_crop,label_canvas,'ScaleFactors',[0.65,0.65,5])
%%
%C(15)
%C(5)
%C(19)
outperm(1025:1196)
%% Visualize any cell array of lists of indexs (local)
tile_list_cell = {outperm(1025:1196),outperm(1218:1301),outperm(1509:1576),outperm(1777:1824),...
    outperm(2007:2040),outperm(2142:2155),outperm(2160:2162),outperm(2164:2222),outperm(2226:2227),outperm(2233:2241)};
thresh = 0.04; 
label_canvas = zeros(400,198,40,'int16');
vol_canvas = zeros(400,198,40,'double');
for i = 1:numel(tile_list_cell)
    cluster_tile_id = int16(tile_list_cell{i});
    for tid = cluster_tile_id
        global_tid = tile_idx(tid); 
        tile_pad = i * int16(tile_space{global_tid} > thresh);
        pos_rng = tile_pos{global_tid};
        label_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
                pos_rng{2}(1):pos_rng{2}(2), ...
                pos_rng{3}(1):pos_rng{3}(2)) = tile_pad;%C_id;
        vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
                pos_rng{2}(1):pos_rng{2}(2), ...
                pos_rng{3}(1):pos_rng{3}(2)) = tile_space{global_tid} + ...
                vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
                pos_rng{2}(1):pos_rng{2}(2), ...
                pos_rng{3}(1):pos_rng{3}(2));
    end
end
%%
volumeViewer(image_vol,label_canvas,'ScaleFactors',[0.65,0.65,5])

%%
for i = 1:numel(tile_list_cell)
fprintf(num2str(tile_list_cell{i}));fprintf('\n')
end
%%
cluster_tile_id = outperm(1025:1196);
vol_canvas = zeros(400,198,40);
for tid = cluster_tile_id
    global_tid = tile_idx(tid);% tid;
    pos_rng = tile_pos{global_tid};
    vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
            pos_rng{2}(1):pos_rng{2}(2), ...
            pos_rng{3}(1):pos_rng{3}(2)) = tile_space{global_tid} + ...
            vol_canvas(pos_rng{1}(1):pos_rng{1}(2), ...
            pos_rng{2}(1):pos_rng{2}(2), ...
            pos_rng{3}(1):pos_rng{3}(2));
end
figure(24);histogram(vol_canvas(vol_canvas>0)) % show the histogram of volume values
%
rX = (1:size(vol_canvas,1))*0.65;
rY = (1:size(vol_canvas,2))*0.65;
rZ = (1:size(vol_canvas,3))*5;
figure(24);hold on 
% p_full
p = patch(isosurface(rY,rX,rZ,vol_canvas,0.035));
isonormals(rY,rX,rZ,vol_canvas,p)
p.FaceColor = C_color_list{i};%'magenta';
p.EdgeColor = 'none';
p.FaceAlpha = 0.3;

p_full = patch(isosurface(rY,rX,rZ,image_vol,0.0035));
isonormals(rY,rX,rZ,image_vol, p_full)
p_full.FaceColor = 'magenta';%C_color_list{i};%'magenta';
p_full.EdgeColor = 'none';
p_full.FaceAlpha = 0.08;

xlabel("X (rost-caudal)")
ylabel("Y (med-lat R-L)")
zlabel("Z (dors-vent)")
set(gca,'Zdir','reverse')
daspect([1 1 1])
% view(3); 
view([-47, 23])
axis tight equal

