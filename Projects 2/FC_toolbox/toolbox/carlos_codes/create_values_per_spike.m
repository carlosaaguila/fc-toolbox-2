function create_values_per_spike(whichPts)
addpath '/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/tools/'
%addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/tools/'
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
run_times_slice={};
fname = {};
gdf_slice = {};
for im = 2 %only want the CAR montage
    for f = 1:length(out.file) %length of file
        for h = 1:length(out.file(f).run) %length of runtimes
            gdf_slice{f}{h} = (out.file(f).run(h).data.montage(im).spikes); %generates gdf for file(1:2), all runtimes, CAR montage.
            run_times_slice{f}{h} = out.file(f).run(h).run_times;
            fname{f} = pt(p).ieeg.file(f).name;
        end
    end
end

gdf={};
all_spikes = {};

%first concatenate all spikes
all_spikes = [];
for im = 2
for f = 1:length(out.file)
    for h = 1:length(out.file(f).run)
            gdf = out.file(f).run(h).data.montage(im).spikes;

            all_spikes = [all_spikes;gdf,...
                repmat(h,size(gdf,1),1),...
                repmat(f,size(gdf,1),1)];
    end
end
end

%pick out 1000 spikes randomly from all files for a single patient. 
n_sp = 1000;
select_spikes = zeros(1000,7);

randi_pathway = '/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/results/mat_output_v2/randi_lists/randi_%s.mat';
filename_randi = sprintf(randi_pathway,pt_name);

if exist(filename_randi) ~= 0
    select_spikes = load(filename_randi,'-mat');
    select_spikes = select_spikes.select_spikes;
else
    for i= 1:n_sp
        sp = randi(size(all_spikes,1));
        f = all_spikes(sp,4);
        h = all_spikes(sp,3);
        fs = out.file(f).run(h).data.fs;
        run_time_slice = out.file(f).run(h).run_times;
        run_start = out.file(f).run(h).run_times(1);
        sp_index = all_spikes(sp,2);          
        sp_time = (sp_index-1)/fs + run_start;        
        sp_ch = all_spikes(sp,1);
        select_spikes(i,:) = [sp_index, sp_ch, h, f, sp_time, run_time_slice];
    end
    save(filename_randi, 'select_spikes');
end
%generate value, fs, channel structures for each patient

values_all = {};
ch_labels_all = {};
fs_all = {};

% Set the pathways
pathway_values = '/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/results/mat_output_v2/values/values_%s.mat';
filename_values = sprintf(pathway_values,pt_name);

pathway_chlabels = '/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/results/mat_output_v2/chlabels/chlabels_%s.mat';
filename_chlabels = sprintf(pathway_chlabels,pt_name);

pathway_fs = '/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/results/mat_output_v2/fs/fs_%s.mat';
filename_fs = sprintf(pathway_fs,pt_name);

if exist(filename_values) ~= 0
    disp('file exists')
    values_all = load(filename_values,'-mat');
    values_all = values_all.values_all;
    ch_labels_all = load(filename_chlabels,'-mat');
    ch_labels_all = ch_labels_all.ch_labels_all;
    fs_all = load(filename_fs,'-mat');
    fs_all = fs_all.fs_all;
    
    if length(values_all) ~= 1000
        for J = length(values_all):1000
            f_all = select_spikes(:,4);
            sp_time_all = select_spikes(:,5);
            data = download_ieeg_data(fname{f_all(J)},login_name,pwfile,[sp_time_all(J)-2,sp_time_all(J)+2],1); % 1 means get lots of data
            disp(J);
            disp('of');
            disp('1000');
            chLabels = data.chLabels;
            values = data.values;
            fs = data.fs;
            % Cleaned labels
            clean_labels = decompose_labels(chLabels,fname{f_all(J)});
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
            fs_all{J} = fs;
            
            savepoints = [500,1000];
            if J == savepoints(1)
                disp('starting save (halfway)-from existing')
                save(filename_values, 'values_all', '-v7.3')
                save(filename_chlabels, 'ch_labels_all')
                save(filename_fs, 'fs_all')
                disp('saved halfway through patient')
            elseif J == savepoints(2)
                disp('starting final save-from existing')
                save(filename_values, 'values_all', '-v7.3')
                save(filename_chlabels, 'ch_labels_all')
                save(filename_fs, 'fs_all')
                disp('saved patient')
            else 
                disp(' no save yet ')
            end
        end
    else
        disp('file is already done, skipping to next patient')
    end
%if files don't exist we got to make em.
else 
    disp('starting new patient!')
    disp('-----------------------')
    for J = 1:1000
            f_all = select_spikes(:,4);
            sp_time_all = select_spikes(:,5);
            data = download_ieeg_data(fname{f_all(J)},login_name,pwfile,[sp_time_all(J)-2,sp_time_all(J)+2],1); % 1 means get lots of data
            disp(J);
            disp('of');
            disp('1000');
            chLabels = data.chLabels;
            values = data.values;
            fs = data.fs;
            % Cleaned labels
            clean_labels = decompose_labels(chLabels,fname{f_all(J)});
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
            fs_all{J} = fs;
            
            savepoints = [500,1000];
            if J == savepoints(1)
                disp('starting save (halfway)- from new')
                save(filename_values, 'values_all', '-v7.3')
                save(filename_chlabels, 'ch_labels_all')
                save(filename_fs, 'fs_all')
                disp('saved halfway through patient')
            elseif J == savepoints(2)
                disp('starting final save- from new')
                save(filename_values, 'values_all', '-v7.3')
                save(filename_chlabels, 'ch_labels_all')
                save(filename_fs, 'fs_all')
                disp('saved patient')
            else 
                disp(' no save yet ')
            end
        end
end

end