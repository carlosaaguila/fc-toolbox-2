addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/tools'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/do_run'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/eeg_processing'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/connectivity_calcs'
addpath '/gdrive/public/USERS/aguilac/Projects/FC_toolbox/toolbox/misc_tools'


locations = fc_toolbox_locs_agui;


whichPts = [99,1,2,3,4,59,43,23]; %Generate montage plots for spikes for these patients

long_run(whichPts);

which_ver = 1;
overwrite = 1;
plot_example_detections1(whichPts,which_ver,overwrite)
