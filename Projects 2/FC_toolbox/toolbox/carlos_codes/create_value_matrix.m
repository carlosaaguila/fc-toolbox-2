function create_value_matrix(whichPts)
%addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/tools/'
addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/tools/'
locations = fc_toolbox_locs_agui;
results_folder = [locations.main_folder,'results/'];
%allout_file = [results_folder, 'all_out/'];
allout_file = '/gdrive/public/USERS/erinconr/projects/fc_toolbox/results/all_out/';%add location to access the PC files generated from long_run
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.main_folder,'data/'];

%load pt file

pt = load([data_folder,'pt.mat']);
pt = pt.pt;

% getting spike information from patient
for p = whichPts
       pt_name = pt(p).name;
end

out = load([allout_file,sprintf('%s_pc.mat',pt_name)]);
out = out.pc;

% generate gdf 
% generate run_times
% generate fnames if multiple
gdf_slice = {};
run_times_slice={};
fname = {};
for im = 2 %only want the CAR montage
    for f = 1:length(out.file) %length of file
        for h = 1:length(out.file(f).run) %length of runtimes
            gdf_slice{f}{h} = (out.file(f).run(h).data.montage(im).spikes); %generates gdf for file(1:2), all runtimes, CAR montage.
            run_times_slice{f}{h} = out.file(f).run(h).run_times;
            fname{f} = pt(p).ieeg.file(f).name;
        end
    end
end

%TEST a slice
%gdf_slice2 = {};
%gdf_slice2{1} = gdf_slice{1}(1,249:253);
%gdf_slice = gdf_slice2;
%run_times_slice2 = {};
%run_times_slice2{1} = run_times_slice{1}(1,249:253);
%run_times_slice= run_times_slice2;

values_all = {};
ch_labels_all = {};
fs_all = {};

for f = 1:length(out.file)
pathway_split = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/split_%s.mat';
filename_split= sprintf(pathway_split,fname{f});
pathway_split2 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/output7.3/split_%s.mat';
filename_split2= sprintf(pathway_split2,fname{f});

if exist(filename_split) ~= 0
    disp('split already exists');
else
pathway1 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/values/values_%s.mat'; %too large
filename1 = sprintf(pathway1,fname{f});
pathway2 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/chlabels/chLabels_%s.mat'; %saved
filename2 = sprintf(pathway2,fname{f});

if exist(filename1) ~= 0
values_all = load(filename1,'-mat');
values_all = values_all.values_all;
%values_all = values_all.values_all;

ch_labels_all = load(filename2,'-mat');
ch_labels_all = ch_labels_all.ch_labels_all;
%ch_labels_all = ch_labels_all.ch_labels_all;

%incorporate a last run type of deal, so that if program crashes due to
%non-server errors, you can pick it up again where you left off

if length(values_all) ~= length(gdf_slice{f})
    for J = length(values_all):length(gdf_slice{f})
        disp('file exists, we will pick up from last save');
        data = download_ieeg_data(fname{f},login_name,pwfile,run_times_slice{f}{J},1); % 1 means get lots of data
        disp(fname{f});
        disp(J);
        disp('of');
        disp(length(gdf_slice{f}));
        chLabels = data.chLabels;
        values = data.values;
        fs = data.fs;
        % Cleaned labels
        clean_labels = decompose_labels(chLabels,fname{f});
        % Find non-intracranial chs
        non_intracranial = find_non_intracranial(clean_labels);
        which_chs = find(~non_intracranial); % channels to do analysis on
        % Reject bad channels
        [bad,~] = identify_bad_chs(values,which_chs,chLabels,fs);
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
        ch_labels_all{J} = curr_labels;
        values_all{J} = run_values;
        %%%%%%%%%%%%
        range = length(;
        for range = range      
            if J == range
                pathway1 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/values/values_%s.mat'; %too large
                pathway2 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/chlabels/chLabels_%s.mat'; %saved
                filename1 = sprintf(pathway1,fname{f});
                filename2 = sprintf(pathway2,fname{f});
                save(filename1, 'values_all','-v7.3');
                save(filename2, 'ch_labels_all');
                disp('saved values and chlabels for this runtime');
            else
                disp('no save yet')
            end
        end
    end
else
    data = download_ieeg_data(fname{f},login_name,pwfile,run_times_slice{f}{1},1);
    disp(fname{f});
    disp('file exists... we load it in');
    fs = data.fs;
end
%

else
    
for d = 1:length(gdf_slice{f})
    
data = download_ieeg_data(fname{f},login_name,pwfile,run_times_slice{f}{d},1); % 1 means get lots of data
disp(fname{f});
disp(d);
disp('of');
disp(length(gdf_slice{f}));
chLabels = data.chLabels;
values = data.values;
fs = data.fs;
fs_all{d} = data.fs;
% Cleaned labels
clean_labels = decompose_labels(chLabels,fname{f});
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
range = floor(linspace(1,length(gdf_slice{f}),30));
    for range = range      
        if d == range
            pathway1 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/values/values_%s.mat'; %too large
            pathway2 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/chlabels/chLabels_%s.mat'; %saved
            filename1 = sprintf(pathway1,fname{f});
            filename2 = sprintf(pathway2,fname{f});
            save(filename1, 'values_all','-v7.3');
            save(filename2, 'ch_labels_all');
            disp('saved values and chlabels for this runtime');
        else 
            disp('no save yet')
        end
    end
end


sprintf('finished %s',fname{f})
pathway1 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/values/values_%s.mat'; %too large
pathway2 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/chlabels/chLabels_%s.mat'; %saved
filename1 = sprintf(pathway1,fname{f});
filename2 = sprintf(pathway2,fname{f});
save(filename1, 'values_all','-v7.3');
save(filename2, 'ch_labels_all');
disp('-----')
disp('FINAL save for values and chlabels for pt');
end

seqs_all = {};
rl = {};
global_coi ={};
leaders = {};

fs = out.file(f).run(1).data.fs;
for k = 1:length(gdf_slice{f})
    [~,rl{k},global_coi{k},~,leaders{k},seqs_all{k}] = build_sequences(gdf_slice{f}{k},length(ch_labels_all{1}),fs); %generates a sequence file for a given gdf
end

% split the values into index ranges, cutting down on the amount of stored
%%values.

%reconstruct gdf from sequences
seqs_concat = {};

for i = 1:length(seqs_all) 
    if isempty(seqs_all{i}) == 1
        disp('...empty sequence cell, proceeding to next...')
    else
        initial = seqs_all{i}{1};
        for k = 2:length(seqs_all{i})
        initial = cat(1,initial,seqs_all{i}{1,k});
        end
    seqs_concat{i} = initial;
    end
end

disp('sequences concatenated')
%parse values above and below each spike in a sequence
%want 1 second around each spike - so 256 each way from the spke peak
%should do it basically 0.5s each way.
value_per_spikeinseq={};

for k = 1:length(seqs_concat)
    for i = 1:length(seqs_concat{k})
        if seqs_concat{k}(i,2)-256 <= 0 %edge condition for the beginning of the values
            ini = values_all{k}(1:seqs_concat{k}(i,2)+256, seqs_concat{k}(i,1))';
            mean1 = nanmean(ini);
            onesbymean = ones(1,513).*mean1;
            onesbymean(:,1:length(ini)) = ini;
            value_per_spikeinseq{k}(i,:) = onesbymean;
            
        elseif seqs_concat{k}(i,2)+256 > (fs*60) %edge condition for the end of the values
            endl = values_all{k}(seqs_concat{k}(i,2)-256:end, seqs_concat{k}(i,1))';
            mean2 = nanmean(endl);
            mean2byones = ones(1,513).*mean2;
            mean2byones(:,end-length(endl)+1:end) = endl;
            value_per_spikeinseq{k}(i,:) = mean2byones;
            
        else
            value_per_spikeinseq{k}(i,:) = values_all{k}(seqs_concat{k}(i,2)-256:seqs_concat{k}(i,2)+256, seqs_concat{k}(i,1))';
        end
    end
end

% create smaller MAT file for python import
split.sequences = seqs_all;
split.values = value_per_spikeinseq;
split.chLabels = ch_labels_all;
split.leaders = leaders;
split.global_coi = global_coi;

pathway_split = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/split_%s.mat';
filename_split= sprintf(pathway_split,fname{f});
pathway_split2 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/output7.3/split_%s.mat';
filename_split2= sprintf(pathway_split2,fname{f});
save(filename_split,'-struct','split');
save(filename_split2,'-struct','split','-v7.3');
disp('saved split into mat_output for pt');

end
end
end