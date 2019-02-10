%% Matt Bachmeier
% ECE 203
% 4/23/2017
% Radio Havana Part 2

%% 4.4.1 Using the DFT

load shortwave.mat
x = raw(:,1) + li*raw(:,2);
X = fftshift(fft(x));
N = length(X);
f_LO = 6000e3;
f_i = 6000e6; % from previous lab

a_raw = x.*exp(-li*2*pi*(1:N)'*(f_i - f_LO)/Fs);
a_raw_X = fftshift(fft(a_raw));
decibels = 20*log10(abs(a_raw_X));
fhat = linspace(-1/2, 1/2, N);
plot(fhat, decibels);
xlabel('Digital frequency');
ylabel('Decibels');
title('Spectrum Before Filtering');

%Setting estimated bandwidth for the station
bandwidth = 4000;
%new filtered signal
new_X = a_raw_X.*(abs(fhat') <= bandwidth/Fs);
decibels = 20*log10(abs(new_X) + 1e-6);

%plot the new filtered signal
plot(fhat, decibels);
xlabel('Digital frequency');
ylabel('Decibels');
title('Spectrum after filtering with DFT');

soundsc(abs(ifft(new_X)), Fs);