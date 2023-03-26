function spike_rate(whichPt)
    addpath '/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/tools/'
    %addpath '/Users/carlosaguila/Desktop/SSH MODEL/Projects/FC_toolbox/tools/'
    locations = fc_toolbox_locs_agui;
    results_folder = [locations.main_folder,'results/'];
    %allout_file = [results_folder, 'all_out/'];
    allout_file = '/mnt/leif/littlab/users/erinconr/projects/fc_toolbox/results/all_out/';%add location to access the PC files generated from long_run
    ieeg_folder = locations.ieeg_folder;
    addpath(genpath(ieeg_folder));
    pwfile = locations.ieeg_pw_file;
    login_name = locations.ieeg_login;
    addpath(genpath(locations.script_folder));
    data_folder = [locations.main_folder,'data/'];
    %{
    gameplan:
    
    Parse through pc.file.run.data.montage.spikes for each patient.
    Tally all the spikes that occur in that run time for each channel.
    Come up with a running tally, then try to get the total time that was 
    looked at. Get Spike Rate per channel as a function of #Spikes/Hour
    
    Now each spike run is 1 minute. Some have none. Some have alot. When we
    look you got to take into account all the runs that were parsed
    through. 
    %}

    pt = load([data_folder,'pt.mat']);
    pt = pt.pt;

    % getting spike information from patient
    for p = whichPt
           pt_name = pt(p).name;
    end
    
    disp(pt_name)
    out = load([allout_file,sprintf('%s_pc.mat',pt_name)]);
    out = out.pc;
    
    %set up an environment
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
    
    spike_rate_matrix = {};
    spike_list = {};
    for f = 1:length(out.file)
        big = 1;
        for h = 1:length(out.file(f).run)
            uniq_chs = [];
            if isempty(gdf_slice{f}{h})
                uniq_chs = [0]; %make a pplaceholder here for those that are empty.
            else
                uniq_chs = unique(gdf_slice{f}{h}(:,1));
            end
            
            yo = 1;
            for i = 1:length(uniq_chs)
                
                if uniq_chs == [0]
                    spike_rate_matrix{f}{h} = {};
                else
                    spike_rate_matrix{f}{h}(yo,1) = sum(gdf_slice{f}{h}(:,1) == uniq_chs(i)); % think to add a new index after f,h called {X} that puts each count/for each uniqch into it's own bracket to correspond with gdf_slice\
                    spike_rate_matrix{f}{h}(yo,2) = uniq_chs(i);
                    spike_list{f}(big,1) = sum(gdf_slice{f}{h}(:,1) == uniq_chs(i));
                    spike_list{f}(big,2) = uniq_chs(i);
                    yo = yo+1;
                    big = big+1;
                end
                
            end
        end
    end
    
    spike_uniq_count = {};
    for f = 1:length(out.file)
    spike_uniq_count{f}(:,2) = unique(spike_list{f}(:,2));
    yo = 1;
        for i = 1:length(spike_uniq_count{f}(:,2))
            spike_uniq_count{f}(yo,1) = sum(spike_list{f}(:,2) == spike_uniq_count{f}(i,2));
            yo = yo+1;
        end
    end
    
    spikerate_perhour = {};
    for f = 1:length(out.file)
        spikerate_perhour{f}(:,1) = spike_uniq_count{f}(:,1) / (length(gdf_slice{f})/60);
        spikerate_perhour{f}(:,2) = spike_uniq_count{f}(:,2);
    end
    
    
    %SAVE
    spike_rate.unique_counts = spike_uniq_count;
    spike_rate.perhour = spikerate_perhour;
    spike_rate.ALL = spike_rate_matrix;
    
    pathway_values = '/mnt/leif/littlab/users/aguilac/Projects/Spike_Rates/spike_rate_%s.mat';
    filename_values = sprintf(pathway_values,pt_name);
    save(filename_values,'-struct','spike_rate');
    
end

