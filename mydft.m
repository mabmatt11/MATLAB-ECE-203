function X = mydft(x)
N=length(x);
for k=1:N
    X(k) = 0;
    for n = 1:N
        X(k) = X(k) + x(n)*exp(-j*2*pi*(k-1)/N*(n-1));
    end
end