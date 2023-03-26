addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/tools'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/do_run'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/eeg_processing'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/connectivity_calcs'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/misc_tools'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/elec_info'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/ccep_stuff'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/create_pt_struct'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/ieeg_stuff'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/net_corr'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/plotting'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/seizure_info'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/elec_info'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/spike_detector'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/spike_sequence_stuff'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/time_grouping'

locations = fc_toolbox_locs_agui;

whichPts = [2,3,4,59,43,23]; %Generate montage plots for spikes for these patients

long_run(whichPts);

which_ver = 1;
overwrite = 1;
plot_example_detections(whichPts,which_ver,overwrite)

