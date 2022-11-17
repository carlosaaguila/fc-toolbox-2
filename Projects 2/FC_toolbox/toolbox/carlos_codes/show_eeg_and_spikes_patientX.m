%%find data all eeg spikes for a channel 


addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/data'; %where the pc data is found.
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/results/all_out'; %where the pt data is found.
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/ieeg_stuff' %where ieegdownload is located
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/tools'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/do_run'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/eeg_processing'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/connectivity_calcs'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/misc_tools'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/elec_info'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/ccep_stuff'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/create_pt_struct'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/ieeg_stuff'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/net_corr'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/plotting'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/seizure_info'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/elec_info'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/spike_detector'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/spike_sequence_stuff'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/toolbox/time_grouping'
addpath '/Users/carlosaguila/Desktop/SSH Model/tools/ieeg-matlab-1.14.49/IEEGToolbox'

locations = fc_toolbox_locs_agui;

data_folder = [locations.main_folder,'data/'];

%load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

whichPts = [99];

% getting spike information from patient
for p = whichPts
       pt_name = pt(p).name;
end

results_folder = [locations.main_folder,'results/'];
allout_file = [results_folder, 'all_out/'];

out = load([allout_file,sprintf('%s_pc.mat',pt_name)]);
out = out.pc;

% generate gdf from out file specificed by "p" (index of the 
gdf = {};
for im = 1:2
    for f = 1:length(out.file)
        for h = 1:length(out.file(f).run)
            gdf{h} = (out.file(f).run(h).data.montage(im).spikes);
        end
    end
end

run_times={};
for q = 1:length(out.file(f).run)
    run_times{q} = out.file.run(q).run_times;
end


%generate run_times
%run_times ={};
%for t = 1:size(pt(p).ieeg.file(f).run_times,1)
%    run_times{t} = pt(p).ieeg.file(f).run_times(t,:);
%end
%run_times2 = [run_times{1}(1),run_times{end}(end)];

% generate patients name
fname = pt(p).ieeg.file(f).name;

pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% plot eeg_and_spikes for patient X

%values, chlabels, fs pulled from IEEG download
%for d = 1:length(gdf)
data = download_ieeg_data(fname,login_name,pwfile,run_times{100},1); % 1 means get lots of data
chLabels = data.chLabels;
values = data.values;
fs = data.fs;
h=figure;
show_eeg_and_spikes(values,chLabels,gdf,fs)
destination='/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/results/all_spikes_and_eeg/HUP203/fig';
%saveas(h,[destination,num2str(d)],'jpeg');
%close(figure)
%end
