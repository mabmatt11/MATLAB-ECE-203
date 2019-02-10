%%  Matt Bachmeier
%  ECE 203
%  Lab 02: Introduction to Complex Exponentials
%  2/5/2017

%% SECTION 4.1

function [xx, tt] = one_cos(amp, ff, phase, dur)
tt = 0:(2*pi)/(20*ff):dur;		%-- gives 20 samples per period
xx = amp*cos(ff*tt+ phase);

rr = one_cos(95, 200*pi, pi/5, 0.025);
plot(rr)

%% SECTION 4.2.1

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
subplot(4,1,1)
plot(tt, xx)
title('Plot Sinusoid for Part 5a (Matt Bachmeier)');
xlabel('Time (sec)');
ylabel('Ampilitude');


%% SECTION 5a

fk = [0.5, 0.5, 0.5];
Xk = [2, 2*exp(j*pi*(-1.25)), (1-j)];
fs = 75;
dur = 6;
tstart = -0.5;

[x0, t0] = syn_sin(fk, Xk, fs, dur, tstart);

%% SECTION 5b

% From the graph, the amp is.........
% The first peak is at ............
% The next peak is at.............
% The period is ......................
% The phase shift is -2*pi*timeShift*frequency

timeshift = -0.0733;
frequency = 0.5;
phase = -2*pi*timeshift*frequency

% Approximate phase is ........ using
% the data cursor tool

%% SECTION 5c

% (1 - j) can be represented as sqrt(2)*exp(j*atan(-1))

% Complex amplitude in polar form
x1 = 2;
x2 = 2*exp(j*pi*(-1.25));
x3 = sqrt(2)*exp(j*atan(-1));

% In rectangular form
% totx = 2 + 2*cos(pi*(-1.25)) + j*2*sin(pi*(-1.25)
%	+ sqrt(2)*cos(atan(-1)) + j*sqrt(2)*sin(atan(-1))

% Add totx together
% totx = 1.5858 + j*0.4142

% Back to polar
% totx = 1.6390*exp(j*0.2555)
% amp is 1.6390 and phase is 0.2555 radians

% With MATLAB
totx = x1 + x2 + x3;
abs(totx)
angle(totx)

%% SECTION 6a

xv = 1; %From part c
dt = 1500;

t1(xv) = sqrt((xv)^2 + (dt)^2)/(3*10^8)

%% SECTION 6b

dyr = 900;
dxr = 100;

t2(xv) = (sqrt((dxr)^2 + (dyr - dt)^2) + sqrt((xv - dxr)^2 + (dyr)^2))/(3*10^8)

%% SECTION 6c

% Defining variables
xv = 1; % Functions can't start at 0
dxr = 100;
dyr = 900;
dt = 1500;
f0 = 150*10^6;
T0 = 1/f0;
tt = -2*T0:1*10^-10:2*T0;

% Equations for comparisons
ss1 = cos(2*pi*f0*(tt - t1));
ss2 = cos(2*pi*f0*(tt - t2));
rv = ss1 - ss2;

% Plotting function
subplot(4,1,2), plot(tt, rv)
title('Plot of rv Part 6c (Matt Bachmeier)');
xlabel('Time (sec)');
ylabel('Amp');
amp = max(rv)

%% SECTION 6d

complexRV = exp(j*2*pi*f0*-t1) - exp(j*2*pi*f0*-t2);
complexAmpRV = abs(complexRV);
% MATLAB can handle complex exponents and their addition so
% converting the waves to complex phasors and applying
% their addition will have MATLAB find the magnitude

subplot(4,1,3), plot(tt, complexRV*exp(j*2*pi*f0*tt))
title('Plot of RV with Complex Exponentials (Matt Bachmeier)');
xlabel('Time (sec)');
ylabel('Amp');
% This graph should be the same as the previous

%% SECTION 6e

xv = 0:0.0001:300;
% all other variables remain the same

% next variables are in seconds
allt1 = sqrt((xv).^2 + (dt)^2)/(3*10^8);
allt2 = (sqrt((dxr)^2 + (dyr - dt)^2) + sqrt((xv - dxr).^2 + (dyr)^2))/(3*10^8);

% add all complex amplitudes along the road and find them at each point
allRV = exp(j*2*pi*f0*-allt1) - exp(j*2*pi*f0*-allt2);
allRVamp = abs(allRV);

%% SECTION 6f

% plot the position and peak signals on a graph
subplot(4,1,4), plot(xv, allRVamp)
title('Plot Signal Strength against position (Matt Bcahmeier)');
xlabel('XV (m)');
ylabel('Signal Strength');

%% SECTION 6g

%Our results show that at certain locations there is signal interference to both increase the signal and destroy the signal entirely. The locations where they interfere constructively are the peaks in the fourth plot. The completely destructive interference occurs at the minima on the fourth plot. These instances of interference occur because of the change in path length and therefore time. When the difference in length creates constructive interference the two waves are in phase with each other. Where as when they totally destruct each other they are half a phase out of sync. The further away the car goes from the signal transmitter the less often each peak is reached because their shift is less noticable.