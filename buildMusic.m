function waveform = buildMusic(Nsamp,noteNums,durs,seconds_per_pulse,fs,env,harm)
waveform = zeros(1,Nsamp);
for pulse=1:numel(durs)
    d = durs(pulse);
    if d>0
        k1 = ceil( (pulse -1)*seconds_per_pulse*fs );
        k2 = floor( ((pulse+d)-1)*seconds_per_pulse*fs );
        tt = (k1:k2)/fs;
        xx = key2noteModified(2,noteNums(pulse),tt,env,harm,seconds_per_pulse);
        waveform(k1+1:k2+1) = xx;
    end
end