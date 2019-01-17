function [] = LoadSpikeData(folder)
% LoadSpikeData.m
%  Gavornik Lab open-ephys setup
% go into a folder and extract open ephys data, convert to MATLAB format,
%  get spiking data (for folders with *.spikes files)
%INPUTS:
%        folder - directory to go into, defaults to current directory
%OUTPUTS:
%        a file named CompiledSpikeData_foldername.mat 

if nargin<1
    folder = pwd;
end

cd(folder);

rawFiles = dir('SE*.spikes');
numChans = length(rawFiles);

allSpikeData = cell(numChans,2);

for ii=1:numChans
   [data,timestamps,~] = load_open_ephys_data_faster(rawFiles(ii).name); 
   allSpikeData{ii,1} = data;
   allSpikeData{ii,2} = timestamps;
end

[events,eventTimes,eventInfo] = load_open_ephys_data_faster('all_channels.events');

digEvents = [];
digTimes = [];

diffs = [5;diff(eventTimes)];
inds = find(~(diffs<=0.005));

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

Fs = eventInfo.header.sampleRate;

auxFiles = dir('*A*.continuous');
numAUX = length(auxFiles);

if numAUX>0 % >0
    cutoff = 200;lpFs = 1000;dsLPRate = Fs/lpFs;
    n = 2;
    [lowb,lowa] = butter(n,cutoff/(Fs/2));
    
    [temp,timestamps,~] = load_open_ephys_data_faster(auxFiles(1).name);
    timepoints = length(temp);
    lpLen = length(1:dsLPRate:timepoints);
    lowpassTimes = timestamps(1:dsLPRate:timepoints);
    
    auxData = zeros(lpLen,numAUX);
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
filename = sprintf('CompiledSpikeData_%s.mat',temp(index(end)+1:end));
save(filename,'numChans','rawFiles','events','Fs',...
    'eventTimes','eventInfo','allSpikeData','auxData','lpFs','lowpassTimes');
end
