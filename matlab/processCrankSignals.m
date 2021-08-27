function [rTan_,rRad_,lTan_,lRad_] = processCrankSignals(path,filename,crankLength)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Load file
data = load([path '/' filename '.mat']);

%% Process crank force signals
% Set analog signals
analogLabels = categorical(data.(filename).Analog.Labels);
samplingFactor = data.(filename).Analog.SamplingFactor;

rTrq = data.(filename).Analog.Data(analogLabels == 'RTan',:);

if isempty(data.(filename).Analog.Data(analogLabels == 'RRad',:))
    rRad = data.(filename).Analog.Data(analogLabels == 'RRadial',:);
end

lTrq = data.(filename).Analog.Data(analogLabels == 'LTan',:);
lRad = data.(filename).Analog.Data(analogLabels == 'LRad',:);

% Downsample if required
rTrq = rTrq(1:samplingFactor:end);
rRad = rRad(1:samplingFactor:end);
lTrq = lTrq(1:samplingFactor:end);
lRad = lRad(1:samplingFactor:end);

% Offset force signals at 2.5mV baseline
offset = 2.5;
rTrq = rTrq - offset;
rRad = rRad - offset;
lTrq = lTrq - offset;
lRad = lRad - offset;

% Set signal conversion factors
rangeV = 5; %voltage
rangeNm =1000; %torque Axis1 = 1000, Axis2 = 700
rangeN = 5000; %force Axis1 = 5000, Axis2 = 4000

mv2nm = rangeNm / rangeV;
mv2n = rangeN / rangeV;

% Convert signal voltages to SI Units (Newton meters, Newtons, Radians)
rTrq = rTrq * mv2nm;
rRad = rRad * mv2n;
lTrq = lTrq * mv2nm;
lRad = lRad * mv2n;

% Convert crank torque to tangential force
rTan = rTrq / crankLength;
lTan = lTrq / crankLength;

% Lowpass filter force signals
frameRate = data.(filename).FrameRate;
fpass = 12;
[b,a] = butter(2,fpass/(frameRate/2));

rTan_ = filtfilt(b,a,rTan);
rRad_ = filtfilt(b,a,rRad);
lTan_ = filtfilt(b,a,lTan);
lRad_ = filtfilt(b,a,lRad);

end

