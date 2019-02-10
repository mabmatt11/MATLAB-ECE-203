function [xx,tt] = mychirp( f1, f2, dur, fsamp )
%MYCHIRP generate a linear-FM chirp signal
%
% usage: xx = mychirp( f1, f2, dur, fsamp )
%
% f1 = starting frequency
% f2 = ending frequency
% dur = total time duration
% fsamp = sampling frequency (OPTIONAL: default is 11025)
%
% xx = (vector of) samples of the chirp signal
% tt = vector of time instants for t=0 to t=dur
%
%McClellan, Schafer, and Yoder, Signal Processing First, ISBN 0-13-065562-7.
%Prentice Hall, Upper Saddle River, NJ 07458. 
%c 2003 Pearson Education, Inc.3
if( nargin < 4 ) %-- Allow optional input argument
    fsamp = 11025;
end
dt = 1/fsamp;
tt = 0:dt:dur;
slope = (f2 - f1)/dur;
cosf = 2*pi*(100 + f1*tt + (slope/2)*tt.*tt);
xx = real(1 * exp(j*cosf)); 
