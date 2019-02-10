function [xx,tt] = syn_sin(fk, Xk, fs, dur, tstart)

%SYN_SIN Function to synthesize a sum of cosine waves
% usage:
% [xx,tt] = syn_sin(fk, Xk, fs, dur, tstart)
% fk = vector of frequencies
% (these could be negative or positive)
% Xk = vector of complex amplitudes: Amp*eˆ(j*phase)
% fs = the number of samples per second for the time axis
% dur = total time duration of the signal
% tstart = starting time (default is zero, if you make this input optional)
% xx = vector of sinusoidal values
% tt = vector of times, for the time axis
%
% Note: fk and Xk must be the same length.
% Xk(1) corresponds to frequency fk(1),
% Xk(2) corresponds to frequency fk(2), etc.

if nargin<5, tstart = 0, end
tt = tstart:1/fs:dur;
xx = Xk(1)*exp(j*2*pi*fk(1)*tt);

for k=2:length(fk)
	xx = xx + Xk(k)*exp(j*2*pi*fk(k)*tt);
end
if length(Xk) ~= length(fk)
	error('error', 'Xk ~= fk')
end
%subplot(4,1,1)
%plot(tt, xx)
%title('Plot Sinusoid for Part 5a (Matt Bachmeier)');
%xlabel('Time (sec)');
%ylabel('Ampilitude');
end