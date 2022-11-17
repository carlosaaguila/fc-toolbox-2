%%
for d = 1:length(gdf)
data = download_ieeg_data(fname,login_name,pwfile,run_times{100},1); % 1 means get lots of data
chLabels = data.chLabels;
values = data.values;
fs = data.fs;
h=figure;
show_eeg_and_spikes(values,chLabels,gdf,fs)
destination='/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/results/all_spikes_and_eeg/HUP203/fig';
saveas(h,[destination,num2str(d)],'jpeg');
close(figure)
end
%%
start1 = seqs{1}(1,2); %retrieves times for sequence
end1 = seqs{1}(end,2);
data = download_ieeg_data(fname,login_name,pwfile,run_times{100},1);
values = data.values;
fs = data.fs;
chLabels = data.chLabels;
%%
gdf_test1 = detector_alt(values,fs,11);
show_eeg_and_spikes(values,chLabels,gdf_test1,fs);
%creating sequences and then running detector_alt doens't seem to work even
%at lower tmul values. 
%%
