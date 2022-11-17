n_sp=50;
surround=1;

% initialize figure
figure
set(gcf,'position',[0 0 1400 1000])
tiledlayout(ceil(n_sp/5),5,'tilespacing','tight','padding','tight');

% Loop over spikes
for i = 1:n_sp

%% Randomly pick spike
sp = randi(size(gdf_tmul19,1));
sp_ch = gdf_tmul19(sp,1);
sp_index = gdf_tmul19(sp,2);

%% Plot data
nexttile
plot(linspace(0,surround*2,size(cat1,1)),cat1(:,sp_ch),'linewidth',1);
hold on
plot(surround,cat1(round(sp_index),sp_ch),'o','markersize',10)
    title(sprintf('Spike %d %1.1f s %s file %d',...
        sp),'fontsize',10)

end