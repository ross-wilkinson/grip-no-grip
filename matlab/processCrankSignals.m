function [rTan_,rRad_,lTan_,lRad_,angleR,angleL] = processCrankSignals(path,filename,crankLength,angleSide)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Load file
data = load([path '/' filename '.mat']);

%% Process crank force and angle signals
% Set analog signals
analogLabels = categorical(data.(filename).Analog.Labels);

rTrq = data.(filename).Analog.Data(analogLabels == 'RTan',:);
rRad = data.(filename).Analog.Data(analogLabels == 'RRadial',:);
lTrq = data.(filename).Analog.Data(analogLabels == 'LTan',:);
lRad = data.(filename).Analog.Data(analogLabels == 'LRad',:);
angle = data.(filename).Analog.Data(analogLabels == 'Crank',:);

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
rangeR = 2*pi; %angle

mv2nm = rangeNm / rangeV;
mv2n = rangeN / rangeV;
mv2r = rangeR / rangeV;

% Convert signal voltages to SI Units (Newton meters, Newtons, Radians)
rTrq = rTrq * mv2nm;
rRad = rRad * mv2n;
lTrq = lTrq * mv2nm;
lRad = lRad * mv2n;
angle = angle * mv2r;

% Convert crank torque to tangential force
rTan = rTrq / crankLength;
lTan = lTrq / crankLength;

% Switch crank signal to right side if necessary
switch angleSide
    case 'right'
        shift = (360/100) * 8;
        angleR = angle + deg2rad(shift);
        angleR(angleR>(2*pi)) = angleR(angleR>(2*pi)) - (2*pi);
        angleOpp = angleR - pi;
        angleOpp(angleOpp<0) = angleOpp(angleOpp<0) + (2*pi);
        angleL = angleOpp;
    case 'left'
        shift = (360/100) * 8;
        angleL = angle + deg2rad(shift);
        angleL(angleL>(2*pi)) = angleL(angleL>(2*pi)) - (2*pi);
        angleOpp = angleL - pi;
        angleOpp(angleOpp<0) = angleOpp(angleOpp<0) + (2*pi);
        angleR = angleOpp;
end

% Lowpass filter force signals
frameRate = data.(filename).FrameRate;
fpass = 12;
[b,a] = butter(2,fpass/(frameRate/2));

rTan_ = filtfilt(b,a,rTan);
rRad_ = filtfilt(b,a,rRad);
lTan_ = filtfilt(b,a,lTan);
lRad_ = filtfilt(b,a,lRad);

end

