function show_eeg_spatial_color(values,fs,cat_ptx_rank,spikeStart,spike)%colormap)


figure
set(gcf,'position',[10 10 1200 800])
dur = size(values,1)/fs;
nchs = 32;%length(labels);

offset = 0;
ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for ich = 1:32
    plot(linspace(0,dur,size(values,1)),values(:,cat_ptx_rank(ich)) - offset,'Color','k')%,colormap(ich,:));
    hold on
    
    for i=1:length(spikeStart{1})
        for j=1:length(spikeStart{1}{1})
            C = spikeStart{1}{1}{j};
            index = find([C{:}] == spike);
        end
    end
    plot(index/fs,1,'o')
    values(spike-(2*2000))

    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(values(:,cat_ptx_rank(ich)));
    
    %text(dur+0.05,ch_bl(ich),sprintf('%s',labels{ich}),'fontsize',20)
    
    if ich < nchs
        
        offset = offset - (min(values(:,cat_ptx_rank(ich))) - max(values(:,cat_ptx_rank(ich+1))));
        
    end
    
end

xlabel('Time (seconds)')
ylabel('Voltage (\muV)')
set(gca,'fontsize',20)

end