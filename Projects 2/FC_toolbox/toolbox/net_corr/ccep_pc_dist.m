function all = ccep_pc_dist(pc,ccep)

% change wrap method

%% Parameters
% Which montage (1 = bipolar, 2 = car)
m = 1;

% Binarize the ccep array. I think this is reasonable because I get a lot
% of very high ccep responses occasionally that will drown out other real
% responses.
do_bin = 0;

% force matrices to be symmetric (exclude response chs that are not stim
% chs). I think this is reasonable to do because if I don't do this then I
% would expect that response chs closer to stim chs would have a higher
% in-degree, creating a bias if I am trying to generalize about indegree
% for different anatomical locations.
do_sym = 1; 

%% Get file locs
locations = fc_toolbox_locs_agui;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'corrs/combined_nets/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Pt struct
data_folder = [locations.main_folder,'data/'];
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% find corresponding p in pt
name = pc.name;
for p = 1:length(pt)
    if strcmp(name,pt(p).name)
        break
    end
end


%% Get pc info
out = net_over_time(pc);
out = reconcile_files(out);
nruns = 0;
for f = 1:length(pc.file)
    nruns = nruns + length(pc.file(f).run);
end
block_dur_secs = diff(pc.file(1).run(1).block_times);
run_dur_secs = diff(pc.file(1).run(1).run_times);

%% Get all base labels
% pc
base_pc_labels = out.all_labels;

% ccep
ccep_labels = ccep.chLabels;

%% Find intersecting labels and corresponding indices
% we need to remember to restrict all future things based on these indices
[indices,intersecting_set] = find_intersecting_idx({base_pc_labels,ccep_labels});
pc_idx = indices{1};
ccep_idx = indices{2};

%% Get avg pc network
%{
pc_net = out.montage(m).net;


% Take average over all times for net
pc_net = nanmean(pc_net,2);

% unwrap
pc_net= wrap_or_unwrap_adjacency_fc_toolbox(pc_net);
pc_net(logical(eye(size(pc_net)))) = nan;
nan_chs = sum(~isnan(pc_net),1) == 0;
%}

% Get all pcnet
all_pc = wrap_or_unwrap_adjacency_fc_toolbox(out.montage(m).net);
pc_net = nanmean(all_pc,3);
ns_time = squeeze(nansum(all_pc,2));

% restrict indices
pc_labels = out.montage(m).labels(pc_idx);
pc_net = pc_net(pc_idx,pc_idx);
ns_time = ns_time(pc_idx,:);
% remove nans
%{
pc_net(nan_chs,:) = [];
pc_net(:,nan_chs) = [];
pc_labels(nan_chs) = [];
%}

%% Get ccep stuff
% network and labels
ccep_net = ccep.A;
ccep_labels_bipolar = ccep.bipolar_labels;
ccep_labels_car = ccep.chLabels;

% put labels into same format as pc labels
empty_labels = cellfun(@isempty,ccep_labels_bipolar);
ccep_labels_bipolar(empty_labels) = {'-'};
ccep_labels_car = cellfun(@(x) [x,'-CAR'],ccep_labels_car,'UniformOutput',false);

% Identify stim and response chs
stim_chs = ccep.ch_info.stim_chs;
%response_chs = ccep.ch_info.response_chs;

if m == 1
    ccep_labels = ccep_labels_bipolar;
elseif m == 2
    ccep_labels = ccep_labels_car;
end
orig_ccep_labels = ccep_labels;


ccep_labels = ccep_labels(ccep_idx);
ccep_net = ccep_net(ccep_idx,ccep_idx);
stim_chs = stim_chs(ccep_idx);
%{
stim_labels = ccep_labels(stim_chs);
response_labels = ccep_labels(response_chs);

% Restrict matrix
ccep_net = ccep_net(response_chs,stim_chs);

% calculate indegree and outdegree
outdegree = nansum(ccep_net,1); % returns 1 x nchs column vector, 1 for each stim (but need to reduce chs!)
indegree = nansum(ccep_net,2); % returns nchs x 1 row vector
outdegree(sum(~isnan(ccep_net),1)==0) = nan;
indegree(sum(~isnan(ccep_net),2)==0) = nan;

% Make non-stim or non-response chs nans
outdegree(~stim_chs) = nan;
outdegree = outdegree';
indegree(~response_chs) = nan;
%}

%% Get distance stuff
% Reconcile locs and anatomy with ieeg labels
locs = pt(p).elecs(1).locs;
loc_labels = decompose_labels(pt(p).elecs(1).elec_names);
locs = reconcile_locs_ieeg(base_pc_labels,loc_labels,locs,[]);

% locs
[~,~,bipolar_labels,chs_in_bipolar,~,mid_locs,~] =...
    bipolar_montage(nan(1,length(base_pc_labels)),base_pc_labels,1:length(base_pc_labels),locs,[]);
% Made it to here!!!!!!!!!!!

[~,car_labels] = car_montage(nan(1,length(base_pc_labels)),1:length(base_pc_labels),base_pc_labels);

if m ==1
    labels = bipolar_labels;
    out_locs = mid_locs;
else
    out_locs = locs;
    labels = car_labels;
end

% restrict idx
out_locs = out_locs(pc_idx,:);
labels = labels(pc_idx);


% Make inv dist matrix
A = inverse_dist(out_locs);
A = wrap_or_unwrap_adjacency_fc_toolbox(A);


%% Reconcile labels
indices = reconcile_labels({ccep_labels,pc_labels,labels});

% Restrict matrices to matching labels
ccep_labels = ccep_labels(indices{1});
ccep_net = ccep_net(indices{1},indices{1});
stim_chs = stim_chs(indices{1});

pc_labels = pc_labels(indices{2});
pc_net = pc_net(indices{2},indices{2});
ns_time  = ns_time(indices{2},:);

labels = labels(indices{3});
A = A(indices{3},indices{3});


%% Make sure labels match
%
if ~isequal(pc_labels,ccep_labels) || ~isequal(pc_labels,labels)
    error('labels do not match')
end
%}

%% Binarize CCEP
if do_bin
    ccep_net(~isnan(ccep_net)) = 1;
end

%% Force each matrix to be symmetric
if do_sym
    ccep_net = ccep_net(stim_chs,stim_chs);
    A = A(stim_chs,stim_chs);
    pc_net = pc_net(stim_chs,stim_chs);
    ns_time = ns_time(stim_chs,:);
    
    ccep_labels = ccep_labels(stim_chs);
    pc_labels = pc_labels(stim_chs);
    labels = labels(stim_chs);
end

%% Again, make sure labels match
%
if ~isequal(pc_labels,ccep_labels) || ~isequal(pc_labels,labels)
    error('labels do not match')
end

%% Get the networks out
all.net.dist.data = A;
all.net.dist.name = '1/{dist^2}';

all.net.pc.data = pc_net;
all.net.pc.name = 'EEG Pearson correlation';
all.ns_time = ns_time;

all.net.ccep.data = ccep_net;
all.net.ccep.name = 'CCEP';
all.net.ccep.x = 'Stim';
all.net.ccep.y = 'Response';
all.net.ccep.orig_labels = orig_ccep_labels;
if isfield(ccep,'is_soz')
    all.net.ccep.is_soz = ccep.is_soz;
end

all.labels = labels;
all.pt_name = pc.name;
all.block_dur_secs = block_dur_secs;
all.run_dur_secs = run_dur_secs;
all.nruns = nruns;

%% Save file
save([out_folder,pc.name,'.mat'],'all')
%{

%% Do correlations
% Correlations across columns
outdegree = nansum(ccep_net,1);
ns_pc_cols = nansum(pc_net,1);
ns_dist_cols = nansum(A,1);
out_pc_r = corr(outdegree',ns_pc_cols','rows','pairwise','type',corr_type);
out_dist_r = corr(outdegree',ns_dist_cols','rows','pairwise','type',corr_type);
pc_dist_r_cols = corr(ns_dist_cols',ns_pc_cols','rows','pairwise','type',corr_type);

% Correlations across rows
indegree = nansum(ccep_net,2);
ns_pc_rows = nansum(pc_net,2);
ns_dist_rows = nansum(A,2);
in_pc_r = corr(indegree,ns_pc_rows,'rows','pairwise','type',corr_type);
in_dist_r = corr(indegree,ns_dist_rows,'rows','pairwise','type',corr_type);
pc_dist_r_rows = corr(ns_dist_rows,ns_pc_rows,'rows','pairwise','type',corr_type);

if abs(pc_dist_r_rows - pc_dist_r_cols) > 1e-4, error('what'); end
if any(abs(ns_pc_cols'-ns_pc_rows) > 1e-4), error('what'); end

%% NS time correlations
%
ns_time_dist_r = corr(ns_dist_rows,ns_time,'rows','pairwise','type',corr_type);
ns_time_in_r = corr(indegree,ns_time,'rows','pairwise','type',corr_type);
ns_time_out_r = corr(outdegree',ns_time,'rows','pairwise','type',corr_type);




%% Plot
figure
set(gcf,'position',[10 10 1300 1000])
tiledlayout(3,3,'tilespacing','tight','padding','tight')

nexttile
turn_nans_gray(A)
%xticks(1:length(labels))
yticks(1:length(labels))
%xticklabels(labels)
xticklabels([])
yticklabels(labels)

title('Inverse distance squared network')
set(gca,'fontsize',15)

nexttile
turn_nans_gray(pc_net)
%{
xticks(1:length(pc_labels))
yticks(1:length(pc_labels))
xticklabels(pc_labels)
yticklabels(pc_labels)
%}
xticklabels([])
yticklabels([])
title('Pearson correlation network')
set(gca,'fontsize',15)

nexttile
turn_nans_gray(ccep_net)
%{
xticks(1:length(ccep_labels))
yticks(1:length(ccep_labels))
xticklabels(ccep_labels)
yticklabels(ccep_labels)
%}
xticklabels([])
yticklabels([])
xlabel('Stim')
ylabel('Response')
title('CCEP network')
set(gca,'fontsize',15)

nexttile
plot(ns_dist_rows,ns_pc_rows,'o','linewidth',2,'markersize',10)
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',pc_dist_r_rows),'verticalalignment','top',...
    'fontsize',15)
xlabel('Inverse distance node strength')
ylabel('PC node strength')
set(gca,'fontsize',15)

nexttile
plot(ns_pc_rows,indegree,'o','linewidth',2,'markersize',10)
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',in_pc_r),'verticalalignment','top',...
    'fontsize',15)
xlabel('PC node strength')
ylabel('Indegree')
set(gca,'fontsize',15)

nexttile
plot(ns_dist_rows,indegree,'o','linewidth',2,'markersize',10)
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',in_dist_r),'verticalalignment','top',...
    'fontsize',15)
xlabel('Inverse distance node strength')
ylabel('Indegree')
set(gca,'fontsize',15)

nexttile
axis off

nexttile
plot(ns_pc_rows,outdegree,'o','linewidth',2,'markersize',10)
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',out_pc_r),'verticalalignment','top',...
    'fontsize',15)
xlabel('PC node strength')
ylabel('Outdegree')
set(gca,'fontsize',15)

nexttile
plot(ns_dist_rows,outdegree,'o','linewidth',2,'markersize',10)
yl = ylim;
xl = xlim;
text(xl(1),yl(2),sprintf('r = %1.2f\n',out_dist_r),'verticalalignment','top',...
    'fontsize',15)
xlabel('Inverse distance node strength')
ylabel('Outdegree')
set(gca,'fontsize',15)

%
figure
set(gcf,'position',[10 10 1300 800])
tiledlayout(2,3,'tilespacing','tight','padding','tight')

% Raster of ns over time
nexttile([1 3])
turn_nans_gray(ns_time)

nexttile
plot(ns_time_dist_r)
title('NS-dist correlation over time')

nexttile
plot(ns_time_in_r)
title('NS-indegree correlation over time')

nexttile
plot(ns_time_out_r)
title('NS-outdegree correlation over time')
%}

end