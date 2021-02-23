function [] = LoadFullData(folder)
% LoadData.m
%  Gavornik Lab open-ephys setup
% go into a folder and extract open ephys data, convert to MATLAB format,
%  lowpass filter and downsample to get LFP, allow option for bandpass
%  filtering to get spiking data
%INPUTS:
%        folder - directory to go into, defaults to current directory
%OUTPUTS:
%        a file named CompiledRawData_foldername.mat 

if nargin<1
    folder = pwd;
end

cd(folder);

[events,eventTimes,eventInfo] = load_open_ephys_data_faster('all_channels.events');

rawFiles = dir('*CH*.continuous');
chansPerTrode = length(rawFiles);numChans = chansPerTrode;
allData = cell(chansPerTrode,1);

for ii=1:chansPerTrode
   [data,timestamps] = load_open_ephys_data_faster(rawFiles(ii).name);
   allData{ii} = data;
end
clear data;

totalTime = max(timestamps)-min(timestamps)+10/1000;
startTime = min(timestamps);

timepoints = length(timestamps);

try 
    fullData = zeros(timepoints,chansPerTrode);
    
    for ii=1:chansPerTrode
        fullData(:,ii) = allData{ii};
    end
    clear allData;
    
    writemda(fullData',sprintf('raw.full.mda'),'float64');
    clear fullData;
catch
    chansPerTrode = 6;
    
    numIter = ceil(numChans/chansPerTrode);
    inds = 1:chansPerTrode;
    for jj=1:numIter
        fullData = zeros(timepoints,chansPerTrode);
        for ii=1:chansPerTrode
            fullData(:,ii) = allData{inds(ii)};
        end
        
        writemda(fullData',sprintf('raw.full%d.mda',jj),'float64');
        clear fullData;
        if jj<numIter-1
            inds = inds+chansPerTrode-1;
        else
            inds = numChans-chansPerTrode+1:numChans;
        end
    end
    clear allData;
end

Fs = eventInfo.header.sampleRate;
% save('RecordingInfo.mat','totalTime','startTime','Fs','timestamps','chansPerTrode');

digEvents = [];
digTimes = [];

diffs = [5;diff(eventTimes)];
inds = find(~(diffs<=0.01));

for ii=1:length(inds)
    digTimes = [digTimes;eventTimes(inds(ii))];
    
    if ii==length(inds)
        currentEvents = inds(ii):length(eventTimes);
        currentEvents = unique(events(currentEvents));
    else
        currentEvents = inds(ii):inds(ii+1)-1;
        currentEvents = unique(events(currentEvents));
    end
    
    binrep = zeros(1,7);
    
    if sum(currentEvents)==0
    
    else
        binrep(currentEvents) = 1;
    end

    digEvents = [digEvents;bi2de(binrep)];
end

events = digEvents;
eventTimes = digTimes;


% notch, low pass, downsample to get LFP
% bandpass to get spiking
cutoff = 200;lpFs = 1000;dsLPRate = Fs/lpFs;
lpLen = length(1:dsLPRate:timepoints);
n = 2;
[lowb,lowa] = butter(n,cutoff/(Fs/2));

lowpassTimes = timestamps(1:dsLPRate:timepoints);

clear timestamps;

auxFiles = dir('*A*.continuous');
numAUX = length(auxFiles);

if numAUX>0 % >0
    auxData = zeros(lpLen,numAUX);
%     d = designfilt('bandpassiir','FilterOrder',6,'HalfPowerFrequency1',2,...
%     'HalfPowerFrequency2',50,'SampleRate',lpFs);
    for ii=1:numAUX
        [temp,~,~] = load_open_ephys_data_faster(auxFiles(ii).name);
        if ~isempty(regexp(auxFiles(ii).name,'ADC1','once'))
    %         temp = filtfilt(notchb,notcha,temp);
            temp = filtfilt(lowb,lowa,temp);
    %         temp = filtfilt(d,temp);
            auxData(:,ii) = temp(1:dsLPRate:timepoints);
        elseif ~isempty(regexp(auxFiles(ii).name,'ADC2','once'))
            threshold = 0.75;
            inds = find(temp(1:end-1)<threshold & temp(2:end)>threshold);
            newinds = ceil(inds./(Fs/lpFs));
            auxData(newinds,ii) = 1;
        elseif ~isempty(regexp(auxFiles(ii).name,'ADC3','once'))
            threshold = 4;
            inds = find(temp(1:end-1)<threshold & temp(2:end)>threshold);
            newinds = ceil(inds./(Fs/lpFs));
            auxData(newinds,ii) = 1;
        end
    end

    for ii=1:numAUX
        if ~isempty(regexp(auxFiles(ii).name,'ADC1','once'))
           temp = auxData(:,ii);
           temp(end-100:end) = mean(temp(end-200:end-101));
           
           checkTime = 0.5*lpFs;
           tempMov = zeros(lpLen,1);
           
           for jj=1:lpLen
               begin = max(jj-checkTime/2,1);
               finish = min(jj+checkTime/2,lpLen);
               y = fft(temp(begin:finish));
               fftLen = ceil(length(y)/2);
               freqs = linspace(0,lpFs/2,fftLen);
               lowInd = 3;
               [~,highInd] = min(abs(20-freqs));
               y = y(1:fftLen);
               power = log(y.*conj(y));
               tempMov(jj) = mean(power(lowInd:highInd));
           end
           %tempMov = tempMov-smooth(tempMov,30*lpFs);
           try
               gm = fitgmdist(tempMov,3);
               [mu,ind] = min(gm.mu);
               sigma = squeeze(gm.Sigma);
               sigma = sigma(ind);
               baseline = mu+3*sigma;
               moveSignal = max(tempMov-baseline,0);
               %            moveSignal = smooth(moveSignal,200);
               %            moveSignal(moveSignal<0.01) = 0;
               auxData(:,ii) = moveSignal;
           catch
               baseline = mean(tempMov)+std(tempMov);
               moveSignal = max(tempMov-baseline,0);
               auxData(:,ii) = moveSignal;
           end
%            figure;plot(temp(1000:50000));
        end
    end
else
   auxData = [];
end

temp = pwd;
index = regexp(temp,'/');
filename = sprintf('CompiledRawData_%s.mat',temp(index(end)+1:end));
save(filename,'totalTime','startTime','Fs','chansPerTrode','numChans','rawFiles','events','lpFs',...
        'eventTimes','eventInfo','lowpassTimes','auxData','Fs','timepoints','-v7.3');
end