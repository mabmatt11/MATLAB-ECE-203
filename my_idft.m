function X = my_idft(x)
N = length(x);
for n = 1:N
    X(n) = 0;
    for k = 1:N
        X(n) = X(n) + x(k)*exp(j*2*pi*(k-1)/N*(n-1));
    end
    X(n) = (1/N)*X(n);
end
