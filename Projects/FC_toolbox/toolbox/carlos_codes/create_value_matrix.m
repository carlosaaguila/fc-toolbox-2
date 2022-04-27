function create_value_matrix(whichPts)
%load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

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

run_times_slice = run_times; %(random);
gdf_slice = gdf; %(random);

%
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

pathway1 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/values/values_%s.mat';
pathway2 = '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/results/mat_output/chlabels/chLabels_%s.mat';
filename1 = sprintf(pathway1,fname);
filename2 = sprintf(pathway2,fname);
save(filename1, values_all);
save(filename2, ch_labels_all);

seqs_all = {};
rl = {};
global_coi ={};
leaders = {};
for k = 1:length(gdf_slice)
    [~,rl{k},global_coi{k},~,leaders{k},seqs_all{k}] = build_sequences(gdf_slice{k},length(ch_labels_all{k}),fs_all{k}); %generates a sequence file for a given gdf
end

% split the values into index ranges, cutting down on the amount of stored
%%values.

%reconstruct gdf from sequences
seqs_concat = {};

for i = 1:length(seqs_all) %for the random 20 sequences]
    initial = seqs_all{i}{1};
    for k = 2:length(seqs_all{i})
        initial = vertcat(initial,seqs_all{i}{1,k});
    end
    seqs_concat{i} = initial;
end

%parse values above and below each spike in a sequence
%want 1 second around each spike - so 256 each way from the spke peak
%should do it basically 0.5s each way.
value_per_spikeinseq={};

for k = 1:length(seqs_concat)
    for i = 1:length(seqs_concat{k})
        if seqs_concat{k}(i,2)-256 <= 0 %edge condition for the beginning of the values
            ini = values_all{k}(1:seqs_concat{k}(i,2)+256, seqs_concat{k}(i,1))';
            mean1 = mean(ini);
            onesbymean = ones(1,513).*mean1;
            onesbymean(:,1:length(ini)) = ini;
            value_per_spikeinseq{k}(i,:) = onesbymean;
            
        elseif seqs_concat{k}(i,2)+256 >= 512*60 %edge condition for the end of the values
            endl = values_all{k}(seqs_concat{k}(i,2)-256:512*60, seqs_concat{k}(i,1))';
            mean2 = mean(endl);
            mean2byones = ones(1,513).*mean2;
            mean2byones(:,end-length(endl):end)
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
filename_split= sprintf(pathway_split,fname);
save(filename_split,'-struct','split');

end