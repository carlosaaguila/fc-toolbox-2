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

whichPts = (1:119);

for i = whichPts
    disp('started:')
    sprintf('pt number %i',i)
    create_values_per_spike(i)
    clear pc
end

%we want to create a +_2sec surround spike, for each channel.
%then we want to store information about where we got the spike from
% time of spike, the run, the file, and overall spike index

