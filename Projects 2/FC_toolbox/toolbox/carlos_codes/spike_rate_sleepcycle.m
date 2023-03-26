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
    
    start_times = readtable('/mnt/leif/littlab/users/aguilac/Projects/Spike_Rates/start_times.csv','Format','auto');
    
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
    
    %initialize start_time for each file.
    start_pt = start_times(whichPt+1,:);
    all_starts = start_pt(:,["Var3","Var4","Var5","Var6","Var7"]);
    start_cleaned = all_starts(:,~(ismissing(all_starts)));
    size_start = size(start_cleaned);
    start_array = table2array(start_cleaned);
    [~, ~, ~, H, MN, S] = datevec(start_array{1});
    starttime_in_sec = H*3600+MN*60+S;
    
    %find chlabels
    chlabels_dir = '/mnt/leif/littlab/users/aguilac/Projects/FC_toolbox/results/mat_output_v2/chlabels/chlabels_%s.mat';
    %chlabels_dir = '/Users/carlosaguila/Downloads/chlabels_%s.mat';
    if exist(sprintf(chlabels_dir,pt_name))
        chlabels = load(sprintf(chlabels_dir, pt_name));
        chlabels = chlabels.ch_labels_all;
        chlabels_clean = chlabels{end};
        chlabels_clean = erase(chlabels_clean, '-CAR');
    else
        data = download_ieeg_data(fname{1},login_name,pwfile,run_times_slice{1}{1},1); % 1 means get lots of data
        chLabels = data.chLabels;
        values = data.values;
        fs = data.fs;
        % Cleaned labels
        clean_labels = decompose_labels(chLabels,fname{1});
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
        ch_labels = curr_labels;
        chlabels_clean = erase(ch_labels, '-CAR');
    end
    
    run_time_secs = {};
    for f = 1:length(out.file)
    run_time_secs{f} = cellfun(@(x) add_start(starttime_in_sec,x), run_times_slice{f}, 'UniformOutput', false);
    run_time_secs{f} = cellfun(@(x)x(1), run_time_secs{f});
    end
    
    %WAKE MOD(X,24) between hour 08 and 20
    spike_rate_matrix_WAKE= {};
    spike_list_WAKE= {};
    
    spike_list_SLEEP = {};
    spike_rate_matrix_SLEEP= {};

    for f = 1:length(out.file)
        big_wake = 1;
        big_sleep = 1;
        for h = 1:length(out.file(f).run)
            
            if run_time_secs{f}(h) >= 8 && run_time_secs{f}(h) <= 20 %hours of wake time 8am to 8pm
                
                uniq_chs_wake = [];
                if isempty(gdf_slice{f}{h})
                    uniq_chs_wake = [0]; %make a placeholder here for those that are empty.
                else
                    uniq_chs_wake = unique(gdf_slice{f}{h}(:,1));
                end
            
                yo_wake = 1;
                for i = 1:length(uniq_chs_wake)

                    if uniq_chs_wake == [0]
                        spike_rate_matrix_WAKE{f}{h} = {[999]};
                    else
                        spike_rate_matrix_WAKE{f}{h}(yo_wake,1) = sum(gdf_slice{f}{h}(:,1) == uniq_chs_wake(i)); % think to add a new index after f,h called {X} that puts each count/for each uniqch into it's own bracket to correspond with gdf_slice\
                        spike_rate_matrix_WAKE{f}{h}(yo_wake,2) = uniq_chs_wake(i);
                        spike_list_WAKE{f}(big_wake,1) = sum(gdf_slice{f}{h}(:,1) == uniq_chs_wake(i));
                        spike_list_WAKE{f}(big_wake,2) = uniq_chs_wake(i);
                        yo_wake = yo_wake+1;
                        big_wake = big_wake+1;
                    end
                end
                
            else %hours of SLEEP

                uniq_chs_sleep = [];
                if isempty(gdf_slice{f}{h})
                    uniq_chs_sleep = [0]; %make a placeholder here for those that are empty.
                else
                    uniq_chs_sleep = unique(gdf_slice{f}{h}(:,1));
                end
            
                yo_sleep = 1;
                for i = 1:length(uniq_chs_sleep)

                    if uniq_chs_sleep == [0]
                        spike_rate_matrix_SLEEP{f}{h} = {[999]};
                    else
                        spike_rate_matrix_SLEEP{f}{h}(yo_sleep,1) = sum(gdf_slice{f}{h}(:,1) == uniq_chs_sleep(i)); % think to add a new index after f,h called {X} that puts each count/for each uniqch into it's own bracket to correspond with gdf_slice\
                        spike_rate_matrix_SLEEP{f}{h}(yo_sleep,2) = uniq_chs_sleep(i);
                        spike_list_SLEEP{f}(big_sleep,1) = sum(gdf_slice{f}{h}(:,1) == uniq_chs_sleep(i));
                        spike_list_SLEEP{f}(big_sleep,2) = uniq_chs_sleep(i);
                        yo_sleep = yo_sleep+1;
                        big_sleep = big_sleep+1;
                    end
                end                
                
            end                 
        end                 
    end
                
    %COUNT
    
    spike_uniq_count_SLEEP = {};
    spike_uniq_count_WAKE = {};
    
    for f = 1:length(spike_list_WAKE)
    spike_uniq_count_SLEEP{f}(:,2) = unique(spike_list_SLEEP{f}(:,2));
    spike_uniq_count_WAKE{f}(:,2) = unique(spike_list_WAKE{f}(:,2));
    yo_SLEEP = 1;
    yo_WAKE = 1;
        for i = 1:length(spike_uniq_count_WAKE{f}(:,2)) % WAKE
            spike_uniq_count_WAKE{f}(yo_WAKE,1) = sum(spike_list_WAKE{f}(:,2) == spike_uniq_count_WAKE{f}(i,2));
            %spike_uniq_count_WAKE{f}(:,3) = zeros(length(spike_uniq_count_WAKE{f}(:,2)),1);
            %spike_uniq_count_WAKE{f}(yo_WAKE,3) = (chlabels_clean(spike_uniq_count_WAKE{f}(i,2)));
            yo_WAKE = yo_WAKE+1;
        end
        for i = 1:length(spike_uniq_count_SLEEP{f}(:,2)) % SLEEP
            spike_uniq_count_SLEEP{f}(yo_SLEEP,1) = sum(spike_list_SLEEP{f}(:,2) == spike_uniq_count_SLEEP{f}(i,2));
            %spike_uniq_count_SLEEP{f}(:,3) = zeros(length(spike_uniq_count_SLEEP{f}(:,2)),1);
            %spike_uniq_count_SLEEP{f}(yo_SLEEP,3) = (chlabels_clean(spike_uniq_count_SLEEP{f}(i,2)));
            yo_SLEEP = yo_SLEEP+1;
        end
    end 
    
    spike_uniq_count_WAKE_cell = {};
    spike_uniq_count_SLEEP_cell = {};
    for f = 1:length(spike_list_WAKE)
        spike_uniq_count_WAKE_cell{f} = num2cell(spike_uniq_count_WAKE{f});
        spike_uniq_count_SLEEP_cell{f} = num2cell(spike_uniq_count_SLEEP{f});
        yo_WAKE = 1;
        yo_SLEEP = 1;
        for i = 1:length(spike_uniq_count_WAKE_cell{f})
            spike_uniq_count_WAKE_cell{f}{yo_WAKE,3} = (chlabels_clean{spike_uniq_count_WAKE_cell{f}{i,2}});
            yo_WAKE = yo_WAKE+1;
        end
        for i = 1:length(spike_uniq_count_SLEEP_cell{f})
            spike_uniq_count_SLEEP_cell{f}{yo_SLEEP,3} = (chlabels_clean{spike_uniq_count_SLEEP_cell{f}{i,2}});
            yo_SLEEP = yo_SLEEP+1;
        end   
    end
    

    %PER HOUR
    
    spikerate_perhour_SLEEP = {};
    spikerate_perhour_WAKE = {};
    for f = 1:length(spike_list_SLEEP)
        spikerate_perhour_SLEEP{f}(:,1) = spike_uniq_count_SLEEP{f}(:,1) / (sum(~cellfun('isempty', spike_rate_matrix_SLEEP{f}))/60);
        spikerate_perhour_SLEEP{f}(:,2) = spike_uniq_count_SLEEP{f}(:,2);
    end
    for f = 1:length(spike_list_WAKE)
        spikerate_perhour_WAKE{f}(:,1) = spike_uniq_count_WAKE{f}(:,1) / (sum(~cellfun('isempty', spike_rate_matrix_WAKE{f}))/60);
        spikerate_perhour_WAKE{f}(:,2) = spike_uniq_count_WAKE{f}(:,2);
    end
    
    spikerate_perhour_SLEEP_cell = {};
    spikerate_perhour_WAKE_cell = {};
    for f = 1:length(spike_list_WAKE)
        spikerate_perhour_WAKE_cell{f} = num2cell(spikerate_perhour_WAKE{f});
        spikerate_perhour_SLEEP_cell{f} = num2cell(spikerate_perhour_SLEEP{f});
        yo_WAKE = 1;
        yo_SLEEP = 1;
        for i = 1:length(spikerate_perhour_WAKE_cell{f})
            spikerate_perhour_WAKE_cell{f}{yo_WAKE,3} = (chlabels_clean{spikerate_perhour_WAKE_cell{f}{i,2}});
            yo_WAKE = yo_WAKE+1;
        end
        for i = 1:length(spikerate_perhour_SLEEP_cell{f})
            spikerate_perhour_SLEEP_cell{f}{yo_SLEEP,3} = (chlabels_clean{spikerate_perhour_SLEEP_cell{f}{i,2}});
            yo_SLEEP = yo_SLEEP+1;
        end   
    end
    
    %SAVE
    spike_rate_SLEEP.unique_counts = spike_uniq_count_SLEEP_cell;
    spike_rate_SLEEP.perhour = spikerate_perhour_SLEEP_cell;
    spike_rate_SLEEP.ALL = spike_rate_matrix_SLEEP;
    
    pathway_values = '/mnt/leif/littlab/users/aguilac/Projects/Spike_Rates/sleep/sleep_spike_rate_%s.mat';
    filename_values = sprintf(pathway_values,pt_name);
    save(filename_values,'-struct','spike_rate_SLEEP');
    
    spike_rate_WAKE.unique_counts = spike_uniq_count_WAKE_cell;
    spike_rate_WAKE.perhour = spikerate_perhour_WAKE_cell;
    spike_rate_WAKE.ALL = spike_rate_matrix_WAKE;
    
    pathway_values = '/mnt/leif/littlab/users/aguilac/Projects/Spike_Rates/wake/wake_spike_rate_%s.mat';
    filename_values = sprintf(pathway_values,pt_name);
    save(filename_values,'-struct','spike_rate_WAKE');
    
end

