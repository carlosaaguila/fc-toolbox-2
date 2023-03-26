% Aguila, Carlos 
% Parses through spike detector outputs to generate 

addpath('/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/tools/')
locations = fc_toolbox_locs_agui;
results_folder = [locations.main_folder,'results/'];
allout_file = [results_folder, 'all_out/']; %add location to access the PC files generated from long_run
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.main_folder,'data/'];

whichPts=(1:147);

for i = whichPts
    %start - print out where we are, UNLESS it's HUP119 which has issues.
    disp('started:')
    sprintf('pt number %i',i)
    
    if i == 17
        continue
    end
    
    spike_rate_per_pt = spike_rate(i); % change FUNCTION 
    
    
    disp('ended pt:')
    sprintf('pt number %i',i)
    disp('-------------------')

    %{
    gameplan:
    
    Parse through pc.file.run.data.montage.spikes for each patient.
    Tally all the spikes that occur in that run time for each channel.
    Come up with a running tally, then try to get the total time that was 
    looked at. Get Spike Rate per channel as a function of #Spikes/Hour
    
    Now each spike run is 1 minute. Some have none. Some have alot. When we
    look you got to take into account all the runs that were parsed
    through.
    
    Create function that takes an ID value - then spits back the spike rate
    per channel. Then accumulate all these into a mat file.
    %}
    
end