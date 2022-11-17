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

locations = fc_toolbox_locs_agui;

whichPts = [2]; %Generate montage plots for spikes for these patients

long_run(whichPts);

which_ver = 1;
overwrite = 1;
plot_example_detections(whichPts,which_ver,overwrite)

