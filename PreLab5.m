lambda = 0.6;
f0 = 3;
T = 20;
x = @(t)t.^2.*exp(-lambda*t).*cos(2*pi*f0*t);
fs = 2*f0+1;
tt = 0:1/fs:T;
xvec = x(tt);
subplot(2,1,1)
plot(tt,xvec); grid on
subplot(2,1,2)
t = linspace(0,T,1000);
plot(t,x(t)); grid on
Xvec = mydft(xvec);