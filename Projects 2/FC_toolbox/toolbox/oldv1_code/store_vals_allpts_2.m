%SCRIPT TO USE CREATE_VALUE_MATRIX.M FUNCTION TO OBTAIN VALUES AND
%SEQUENCES FOR EACH PATIENT

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

whichPts = (60:85); %skipped 59 and 56

for i = whichPts
    disp('started:')
    sprintf('pt number %i',i)
    create_value_matrix(i)
    clear pc
end
