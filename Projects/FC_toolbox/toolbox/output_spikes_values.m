% want to create a script that pulls epochs of eeg data from  patient
% around the seq times that are specified. 
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/tools/'
locations = fc_toolbox_locs_agui;
results_folder = [locations.main_folder,'results/'];
allout_file = [results_folder, 'all_out/']; %add location to access the PC files generated from long_run
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.main_folder,'data/'];

%load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

whichPts = [99]; % loads HUP203

% getting spike information from patient
for p = whichPts
       pt_name = pt(p).name;
end

out = load([allout_file,sprintf('%s_pc.mat',pt_name)]);
out = out.pc;

% generate gdf from out file specificed by "p" (index of the 
gdf = {};
for im = 2 %only want the CAR montage
    for f = 1:length(out.file) %length of file
        for h = 1:length(out.file(f).run) %length of runtimes
            gdf{h} = (out.file(f).run(h).data.montage(im).spikes); %generates gdf for file(1), all runtimes, CAR montage.
        end
    end
end

%generates run times
run_times={};
for q = 1:length(out.file(f).run)
    run_times{q} = out.file.run(q).run_times;
end

% generate patients name
fname = pt(p).ieeg.file(f).name;

%% get a sample dataset. lets go for 10 run_times

run_times_slice = run_times;
gdf_slice = gdf;

%%
values_all = {};
ch_labels_all = {};
fs_all = {};

for d = 1:length(gdf_slice)
data = download_ieeg_data(fname,login_name,pwfile,run_times_slice{d},1); % 1 means get lots of data
chLabels = data.chLabels;
values = data.values;
fs = data.fs;
fs_all{d} = data.fs;
% Cleaned labels
clean_labels = decompose_labels(chLabels,fname);
% Find non-intracranial chs
non_intracranial = find_non_intracranial(clean_labels);
which_chs = find(~non_intracranial); % channels to do analysis on
% Reject bad channels
[bad,details] = identify_bad_chs(values,which_chs,chLabels,fs);
which_chs(ismember(which_chs,bad)) = []; % reduce channels to do analysis on
% CAR montage
[car_values,car_labels] = car_montage(values,which_chs,clean_labels);
is_run_car = ismember((1:length(clean_labels))',which_chs);
values = car_values;
is_run = is_run_car;
curr_labels = car_labels;
values = notch_filter(values,fs);
values = bandpass_filter(values,fs);
% make non run channels nans
run_values = values;
run_values(:,~is_run) = nan;
skip = find(~is_run);
ch_labels_all{d} = curr_labels;
values_all{d} = run_values;
end

%%
seqs_all = {};
rl = {};
for k = 1:length(gdf_slice)
    [~,rl{k},~,~,~,seqs_all{k}] = build_sequences(gdf_slice{k},length(ch_labels_all{k}),fs_all{k}); %generates a sequence file for a given gdf
end

%% create large MAT file for python import
split_1.sequences = seqs_all;
split_1.values = values_all;
split_1.chLabels = ch_labels_all;
save('/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/split_1.mat','split_1', '-v7.3');
%% plot sequence
%figure(1)
%show_eeg_and_spikes(values_all{1},ch_labels_all{1}, seqs_all{1}{1}, fs_all{1}) %this plots 1st sequence in our first gdf

%this is to show that gdf lines up with values