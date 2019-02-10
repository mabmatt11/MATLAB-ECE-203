function xboost = MEGAbass(x,AmplificationFactor,BassCutoffFreq,fs)

xboost = x;
for k = 0:N-1;
    if (k*fs/N) <= fcutoff
        xboost(k+1) = xboost(k+1)*Ampfactor;
        xboost(N-k) = xboost(N-k)*Ampfactor;
    end
end