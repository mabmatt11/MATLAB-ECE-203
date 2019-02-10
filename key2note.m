function xx = key2note(X, keynum, dur,fs)
% KEY2NOTE Produce a sinusoidal waveform corresponding to a
%       given piano key number
%
%  usage: xx = key2note (X, keynum, dur)
%
%      xx = the output sinusoidal waveform
%       X = complex amplitude for the sinusoid, X = A*exp(j*phi).
%  keynum = the piano keyboard number of the desired note
%     dur = the duration (in seconds) of the output note
%

if nargin < 4
    fs = 11025;
end
tt = 0:(1/fs):dur;
freq = 440*2^((keynum-49)/12);
no_envelope = real(X*exp(j*2*pi*freq*tt) + 0.75*X*exp(j*2*pi*2*freq*tt) + 0.5*X*exp(j*2*pi*3*freq*tt) + 0.25*X*exp(j*2*pi*4*freq*tt));
length_of_tone = length(no_envelope);
attack = linspace(0,1,length_of_tone*4/32);
delay = linspace(1,.8,length_of_tone*2/32);
sustain = linspace(.8,.75,length_of_tone*21/32);
release = linspace(.75,0,length_of_tone*5/32);

envelope = [attack,delay,sustain,release];

difference = length_of_tone - length(envelope);
envelope = [envelope, zeros(1,difference)];
xx = envelope.*no_envelope;

        