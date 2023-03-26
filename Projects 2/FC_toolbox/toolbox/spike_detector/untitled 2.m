% CODE TO MAKE USE OF EVENT DETECTOR WRITEN BY HANK BINK

%according to the notes in binks .m files
%finds the constraints at where the experiment starts and ends

load('/Users/carlosaguila/Downloads/CAT/Exp1Filt.mat')
load('/Users/carlosaguila/Downloads/CAT/allExpThElecData.mat')
load('/Users/carlosaguila/Downloads/CAT/AllExpTimingSec.mat')
load('/Users/carlosaguila/Downloads/CAT/dtds.mat')
load('/Users/carlosaguila/Downloads/CAT/Exp1Unfilt.mat')

featStartSecExp = [350 310 660 160 470];
featEndSecExp = [2550.85 7245 3227 7252 6929];
featStart= {};
featEnd= {};
for exper=1:5
featStart{exper} = find(allExpTimingSec{exper}>featStartSecExp(exper),1);
featEnd{exper} = find(allExpTimingSec{exper}>featEndSecExp(exper),1);
end

cat1 = allEventData';
%%
%now bring the sample size to the featStart and featEnd
%experiment 1
nchns= size(cat1,2);
for exper=1
    cat1_trim = zeros(featEnd{exper}-featStart{exper}+1,nchns);
    for nchns = 1:(nchns)
        cat1_trim(:,nchns) = cat1((featStart{exper}:featEnd{exper}),nchns);
    end
end

%% tmul30
fs = [];
for f = 1:length(dtds)
    fs(f) = 1/dtds(f) * 1000;
end
gdf_tmul30 = detector_alt(cat1_trim,fs(1),30);

%% plot some of the spikes detected using conrad detector
figure
set(gcf,'position',[0 0 1400 1000])
tiledlayout(ceil(nchns/5),5,'tilespacing','tight','padding','tight');
surround=1;
for i = 1:50
sp = randi(size(gdf_tmul30,1));
sp_index = gdf_tmul30(sp,2);
sp_ch = gdf_tmul30(sp,1);

nexttile
plot(linspace(0,2,size(cat1_trim,1)),cat1_trim(:,sp_ch),'linewidth',1); %% find a way to plot the area surrounding the spike
hold on
xlim([0.995,1.005])
plot(surround,cat1_trim(round(sp_index),sp_ch),'x','markersize',5)
yticklabels([])
set(gca,'fontsize',10)
end

%% tmul = 50
 
fs = [];
for f = 1:length(dtds)
    fs(f) = 1/dtds(f) * 1000;
end
gdf_tmul50 = detector_alt(cat1_trim,fs(1),50);

%% plot tmul 50
figure
set(gcf,'position',[0 0 1400 1000])
tiledlayout(ceil(nchns/5),5,'tilespacing','tight','padding','tight');
surround=1;
for i = 1:50
sp = randi(size(gdf_tmul50,1));
sp_index = gdf_tmul50(sp,2);
sp_ch = gdf_tmul50(sp,1);

nexttile
plot(linspace(0,2,size(cat1_trim,1)),cat1_trim(:,sp_ch),'linewidth',2); %% find a way to plot the area surrounding the spike
hold on
xlim([0.97,1.03])
plot(surround,cat1_trim(round(sp_index),sp_ch),'o','markersize',10)
yticklabels([])
set(gca,'fontsize',10)
end

%% hankbinks' detector way for cat1_trim
evStart={};
evStop = {};
spikeStart ={};
spikeStop = {};
numSpikes={};
sdThresh = 130;
for detChanData = 1:nchns
[evStart{detChanData}, evStop{detChanData}, spikeStart{detChanData}, spikeStop{detChanData}, numSpikes{detChanData}] = eventDetector((cat1(:,detChanData))', dtds(1), featStart{1}(1), featEnd{1}(1), sdThresh);
end


%% 

sdThresh=130;
[evStart1, evStop1, spikeStart1, spikeStop1, numSpikes1] = eventDetector((cat1(:,4))', dtds(1), featStart{1}(1), featEnd{1}(1), sdThresh);
%channel 4 is the one with the "highest" amp

%% 

for i=1:size(cat1,1)
    if cat1(i,:) == allExpThElecData{1}
        disp(i)
    end
end

for i=1:size(cat2,1)
    if cat2(i,:) == allExpThElecData{2}
        disp(i)
    end
end

for i=1:size(cat3,1)
    if cat3(i,:) == allExpThElecData{3}
        disp(i)
    end
end

for i=1:size(cat4,1)
    if cat4(i,:) == allExpThElecData{4}
        disp(i)
    end
end

for i=1:size(cat5,1)
    if cat5(i,:) == allExpThElecData{5}
        disp(i)
    end
end