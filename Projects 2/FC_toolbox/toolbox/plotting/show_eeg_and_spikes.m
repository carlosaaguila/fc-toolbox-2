function [EEG123,time_EEG123,timed_index,spikes_offset_idx] = show_eeg_and_spikes(values,chLabels,gdf,fs)
%NOTES:
%LOOKS LIKE THE CODE DOESN'T NECESSARILY TAKE OUT NAN'S THE OFFSETS ARE
%JUST TO GET THEM ALL ON THE SAME PLOT. 
%THE NAN'S REMAIN BUT DON'T CAUSE AN OFFSET.

%figure
%set(gcf,'position',[62 104 1145 701])

offset = 0;
ch_offsets = zeros(size(values,2),1);
ch_bl = zeros(size(values,2),1);
dur = size(values,1)/fs;

%initialize new variables I want saved for analyses
spikes_offset_idx = {};
timed_index = {};
EEG123 = {};
time_EEG123 = {};

for ich = 1:6%:size(values,2)
    
    EEG123{ich} = values(:,ich)-offset;
    time_EEG123{ich} = linspace(0,dur,size(values,1));
    plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'k');
    
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(values(:,ich));
    hold on
    text(dur+0.05,ch_bl(ich),sprintf('%s',chLabels{ich}))
    
    if ich<size(values,2)
        if ~isnan(min(values(:,ich)) - max(values(:,ich+1)))
            offset = offset - (min(values(:,ich)) - max(values(:,ich+1)));
        end
    end
end

for s = 1:size(gdf,1)
    %index = spikes(s,1);
    index = gdf(s,2);
    
    % convert index to time
    time = index/fs;
    
    ch = gdf(s,1);
    offset_sp = ch_offsets(ch);
    
    value_sp = values(round(index),ch);
    
    plot(time,value_sp - offset_sp,'ro') %plot redcircle where spikes are located
    spikes_offset_idx{s}  = value_sp - offset_sp;
    timed_index{s} = time;
    
end

%pause
%close(gcf)

end