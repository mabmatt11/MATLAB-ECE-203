function xx = key2noteModified(X,keynum,tt,env,harm,seconds_per_pulse)
% KEY2NOTE Produce a sinusoidal waveform corresponding to a
% given piano key number
%
% usage: xx = key2note (X, keynum, dur)
%
% xx = the output sinusoidal waveform
% tt = times at which to compute the output sinusoid
% X = complex amplitude for the sinusoid, X = A*exp(j*phi).
% keynum = the piano keyboard number of the desired note
% env = true to use envelope; false otherwise
% harm = true to add 2nd and 3rd harmonics, false otherwise
% will call envelope.m

f0 = 440; % A above middle C
key0 = 49; % is key 49 on the piano
freq = f0*2^((keynum-key0)/12);
% fprintf('freq of key %g is %g Hz\n',keynum,freq)
higherharmonics = 0;
if harm
    higherharmonics = .4*exp(2j*pi*(2*freq)*tt) + .2*exp(2j*pi*(3*freq)*tt);
end
fundamental = exp(2j*pi*freq*tt);
xx = real( X*( fundamental + higherharmonics ) );
if env
    xx = xx.*envelope(tt,seconds_per_pulse);
end