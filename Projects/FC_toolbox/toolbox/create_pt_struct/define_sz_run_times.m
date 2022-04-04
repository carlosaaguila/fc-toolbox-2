function define_sz_run_times

%% Establish parameters
overwrite = 1;
% Decide frequency and duration of sampling
mini_block = 60; % each run should be 60 seconds
peri_ictal_time = 3600*2; % 2 hours before and two hours after the sz

%% Seed random number generator (so that I get the same thing every time if I re-run)
rng(0)

%% Get file locs
locations = fc_toolbox_locs_agui;
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

for p = 1:length(pt)

    if isempty(pt(p).ieeg)
        continue
    end
    
    if isfield(pt(p).ieeg.file(1),'sz_run_times') == 1 && ...
            ~isempty(pt(p).ieeg.file(1).sz_run_times)
        if overwrite == 0
            fprintf('\nskipping %s\n',pt(p).name);
            continue
        end
        
    end

    fprintf('\nDoing %s\n',pt(p).name);
    
    % Loop over files
    nfiles = length(pt(p).ieeg.file);
    for f = 1:nfiles
        
        % Loop over szs
        if ~isfield(pt(p).ieeg.file(f),'sz_times'), continue; end
        
        sz_times = pt(p).ieeg.file(f).sz_times;
        nszs = size(sz_times,1);
        pt(p).ieeg.file(f).sz_run_times = cell(nszs,1);
        
        for s = 1:nszs
            sbegin = sz_times(s,1);
            send = sz_times(s,2);
            
            % Build run times
            sz_run_times = nan(peri_ictal_time*2/mini_block,2);
            for sr = 1:peri_ictal_time*2/mini_block
                if sr <= peri_ictal_time/mini_block
                    sz_run_times(sr,1) = sbegin - peri_ictal_time + mini_block*(sr-1); % first block starts peri_ictal_time bacl
                    sz_run_times(sr,2) = sbegin - peri_ictal_time + mini_block*(sr); % last block ends at sbegin
                else
                    sz_run_times(sr,1) = send - peri_ictal_time + mini_block*(sr-1); % first block starts at send
                    sz_run_times(sr,2) = send - peri_ictal_time + mini_block*(sr); % last block ends at peri_ictal_time forward
                end
            end
            
            % add
            pt(p).ieeg.file(f).sz_run_times{s} = sz_run_times;
        
        end
    end

    % Save
    save([data_folder,'pt.mat'],'pt');
end