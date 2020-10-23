function [emgRms, emgRaw, channelNames] = processEmgSignals(path,filename,emgChannelNo,window,Fs,lowcut,highcut,order)
%UNTITLED13 Summary of this function goes here
%   Inputs:
%           - window: window size in milliseconds (ms)
%           - Fs: sampling frequency in Hz
%           - lowcut: low cut off frequency in Hz
%           - high cut off frequency in Hz
%           - order: nth order for butterworth filter

% Load file
data = load([path '/' filename '.mat']);

emgRaw = data.(filename).Analog.Data(emgChannelNo,:)';

% Setup band-pass filter
Fnyq = Fs/2; % Normalized by Nyquist frequency
n = order/2; % Order. Note bandpass and bandstop filters will be of order 2n
Wn = [lowcut highcut]/Fnyq; % Create bandpass cutoff vector        
[b,a] = butter(n,Wn); % Creates 2nd order Butterworth filter.

% Remove DC offset
emgOffset = emgRaw - mean(emgRaw);

% Band-pass filter emg signals
emgFilt = filtfilt(b,a,emgOffset); % Bandpass filter signals

% RMS with specified window size
nSamples = window / 1000 * Fs; % convert to seconds then multiply by sampling freq.
movrmsWin = dsp.MovingRMS(nSamples); % 
emgRms = movrmsWin(emgFilt);

analogLabels = categorical(data.(filename).Analog.Labels);
channelNames = cellstr(analogLabels);

end

