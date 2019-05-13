%% Calculate the center of mass
tic
CoM_array = zeros(length(tile_idx),3);
parfor tid=1:length(tile_idx)
    disp(tid)
    global_tid = tile_idx(tid); 
    tile_patch = tile_space{global_tid};
    pos_rng = tile_pos{global_tid};
    [X,Y,Z] = meshgrid(pos_rng{2}(1):pos_rng{2}(2), ...
                pos_rng{1}(1):pos_rng{1}(2),...
                pos_rng{3}(1):pos_rng{3}(2));
     X = single(X); Y = single(Y); Z = single(Z);
     CoM_tile = [sum(X.*tile_patch, 'all'), ...
                sum(Y.*tile_patch, 'all'), ...
                sum(Z.*tile_patch, 'all') ] ./ sum(tile_patch, 'all');
     CoM_array(tid,:) = CoM_tile;
end
toc
save('tile_centerofmass.mat', CoM_array)
%%
tic
Mass_array = zeros(length(tile_idx),1);
parfor tid=1:length(tile_idx)
    global_tid = tile_idx(tid); 
    tile_patch = tile_space{global_tid};
     Mass_array(tid,:) = sum(tile_patch, 'all');
end
toc
%%
ind_array = outperm(1025:1180);
figure(18)
scatter3(0.65 * CoM_array(ind_array, 1), 0.65 * CoM_array(ind_array, 2), 5 * CoM_array(ind_array, 3))
axis equal
%% Cluster Activity Center of Mass Trajectory ! 
% Calculate the center of mass of activity (raw florescent to avoid negative)
act_CoM_traj = CoM_array(ind_array, :)' * raw_traces(tile_idx(ind_array), :) ./ sum(raw_traces(tile_idx(ind_array), :), 1); %zscore_arr(ind_array, :) ./ sum(zscore_arr(ind_array, :), 1);
% CoM_traj = CoM_array(ind_array, :)' * zscore_arr(ind_array, :) ./ sum(zscore_arr(ind_array, :), 1); %zscore_arr(ind_array, :) ./ sum(zscore_arr(ind_array, :), 1);
time_slice = 7500:8000; % around an activity peak! 
figure(18);clf;hold on
scatter3(0.65 * CoM_array(ind_array, 1), 0.65 * CoM_array(ind_array, 2), 5 * CoM_array(ind_array, 3), 'MarkerFaceAlpha', 0.7)
axis equal
comet3(0.65*act_CoM_traj(1, time_slice), 0.65*act_CoM_traj(2, time_slice), 5*act_CoM_traj(3, time_slice))
%% Visualizing the calcium activity on the 3d space

figure(20)
for t_step = time_slice
    scatter3(0.65 * CoM_array(ind_array, 1), ...
        0.65 * CoM_array(ind_array, 2), ...
        5 * CoM_array(ind_array, 3), 50, ...
        zscore_arr(ind_array, t_step), 'filled')
%         raw_traces(tile_idx(ind_array), t_step), 'filled')
    axis equal
    pause(0.1)
end

%% Fancy version 
samp_rate = 20;
do_record =1 ;
cluster_be = [1776, 1847]; 
ind_array =  outperm(cluster_be(1):cluster_be(2));%1782:1884;%1:2271;
time_slice = 6000:12000;%5479:5800;%4164:8620;%6056:19660;%1:24000;
activity_arr = zscore_arr; % raw_traces(tile_idx, :)

cmax = 5;%max(activity_arr(ind_array, time_slice), [], 'all');
cmin = min(activity_arr(ind_array, time_slice), [], 'all');

if do_record
    clear F
    F(length(time_slice)) = struct('cdata',[],'colormap',[]);
    frame_i = 1;
end
figure(22);clf;hold on
h_str = scatter3(0.65 * CoM_array(:, 1), ...
        0.65 * CoM_array(:, 2), ...
        5 * CoM_array(:, 3), 9, 'black','filled', 'MarkerFaceAlpha', 0.3);
h_act = scatter3(0.65 * CoM_array(ind_array, 1), ...
    0.65 * CoM_array(ind_array, 2), ...
    5 * CoM_array(ind_array, 3), 2*Mass_array(ind_array), ...
    activity_arr(ind_array, 1), 'filled', ...
    'MarkerFaceAlpha', 0.9);
axis equal tight
view([-55, 21])
xlabel("X (rost-caud)")
ylabel("Y (med-lat R-L)")
zlabel("Z (dors-vent)")
caxis([cmin, cmax])
ch = colorbar();
set(get(ch,'Label'), 'string', 'zscore');
ax = gca;
ax.NextPlot = 'replaceChildren';
for t_step = time_slice%time_slice
    activity_vec = activity_arr(ind_array, t_step);
    h_act.CDataSource = 'activity_vec'  ; % only update data1!
    refreshdata
    title(sprintf("%.2f s (frame %d)", t_step/samp_rate, t_step))
    drawnow  
    pause(0.05)
    if do_record
        F(frame_i) = getframe(gcf);
        frame_i = frame_i+1;
    end
end
filename = sprintf('Cluster%d-%d_dyn_t%d-%d',cluster_be(1),cluster_be(2),...
    time_slice(1),time_slice(end));
%%
figure()
h = scatter3(0.65 * CoM_array(:, 1), ...
        0.65 * CoM_array(:, 2), ...
        5 * CoM_array(:, 3), 2*Mass_array, 'b', 'MarkerFaceAlpha', 0.6)
%% Write to AVI
for n = 1:length(F)
    frame = F(n);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if n == 1
      imwrite(imind,cm,[filename, '.gif'],'gif', 'Loopcount',inf);
    else
      imwrite(imind,cm,[filename, '.gif'],'gif','WriteMode','append');
    end
end
%% Write
writerObj = VideoWriter([filename,'.avi']);
writerObj.FrameRate = 10;
 % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
%%
figure(22)
movie(F,2)