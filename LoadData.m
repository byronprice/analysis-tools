function [] = LoadData(folder,keepBandpass)
% LoadData.m
%  Gavornik Lab open-ephys setup
% go into a folder and extract open ephys data, convert to MATLAB format,
%  lowpass filter and downsample to get LFP, allow option for bandpass
%  filtering to get spiking data
%INPUTS:
%        folder - directory to go into, defaults to current directory
%        keepBandpass - logical (1 or 0) to save bandpass data or not
%                 if you just want the LFP, then you type 0, defaults to 0
%OUTPUTS:
%        a file named CompiledData_foldername.mat 

if nargin<1
    folder = pwd;
    keepBandpass = 0;
elseif nargin<2
    keepBandpass = 0;
end

cd(folder);

rawFiles = dir('*CH*.continuous');
numChans = length(rawFiles);
[data1,timestamps,~] = loadAndCorrectPhase(rawFiles(1).name,1);

allData = zeros(length(timestamps),numChans);
allData(:,1) = data1;

for ii=2:numChans
   [data,~,~] = loadAndCorrectPhase(rawFiles(ii).name,1); 
   allData(:,ii) = data;
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

[timepoints,numChans] = size(allData);
Fs = eventInfo.header.sampleRate;

% notch, low pass, downsample to get LFP
% bandpass to get spiking
cutoff = 200;lpFs = 1000;dsLPRate = Fs/lpFs;
lpLen = length(1:dsLPRate:timepoints);
bpFs = 10000;dsBPRate = Fs/bpFs;
bpLen = length(1:dsBPRate:timepoints);
n = 2;
[lowb,lowa] = butter(n,cutoff/(Fs/2));

% wo = 60/(Fs/2);
% bw = wo/n;
% [notchb,notcha] = iirnotch(wo,bw);

d = designfilt('bandpassiir','FilterOrder',n,'HalfPowerFrequency1',cutoff+100,...
    'HalfPowerFrequency2',5000,'SampleRate',Fs);

bandpassData = zeros(bpLen,numChans);
lowpassData = zeros(lpLen,numChans);
lowpassTimes = timestamps(1:dsLPRate:timepoints);
bandpassTimes = timestamps(1:dsBPRate:timepoints);
for ii=1:numChans
%     temp = filtfilt(notchb,notcha,allData(:,ii));
    temp = filtfilt(lowb,lowa,allData(:,ii));
    lowpassData(:,ii) = temp(1:dsLPRate:timepoints);
    
    temp = filtfilt(d,allData(:,ii));
    bandpassData(:,ii) = temp(1:dsBPRate:timepoints);
end

auxFiles = dir('*A*.continuous');
numAUX = length(auxFiles);

if numAUX>0
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
            newinds = round(inds./(Fs/lpFs));
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
           gm = fitgmdist(tempMov,3);
           [mu,ind] = min(gm.mu);
           sigma = squeeze(gm.Sigma);
           sigma = sigma(ind);
           baseline = mu+3*sigma;
           moveSignal = max(tempMov-baseline,0);
%            moveSignal = smooth(moveSignal,200);
%            moveSignal(moveSignal<0.01) = 0;
           auxData(:,ii) = moveSignal;
%            figure;plot(temp(1000:50000));
        end
    end
else
   auxData = [];
end

temp = pwd;
index = regexp(temp,'/');
filename = sprintf('CompiledData_%s.mat',temp(index(end)+1:end));
if keepBandpass == 0
    save(filename,'numChans','rawFiles','events','lpFs',...
        'eventTimes','eventInfo','lowpassTimes','lowpassData','auxData');
else
    save(filename,'numChans','rawFiles','events','eventTimes','eventInfo',...
        'lowpassTimes','lowpassData','bandpassTimes','bandpassData','auxData',...
        'lpFs','bpFs');
end

end
