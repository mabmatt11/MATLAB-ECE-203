%% Matt Bachmeier
% ECE 203
% 4/23/2017
% Radio Havana Part 2

%% 4.4.1 Using the DFT

load shortwave.mat
x = raw(:,1) + 1i*raw(:,2);
X = fftshift(fft(x));
N = length(X);
f_LO = 6000e3;
f_i = 6e6; % from previous lab

a_raw = x.*exp(-1i*2*pi*(1:N)'*(f_i - f_LO)/Fs);
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

%soundsc(abs(ifft(new_X)), Fs);
%The station is much more clear than the end of part 1 of this lab, and
%there is no high pitched noise in the background anymore.

%% 4.4.2(a) Using FIR filters (truncated sinc)

L = 200; %length of the filter
bandwidth = 4000;
fm = bandwidth/Fs;

%Make the impulse response of truncated sinc
k = 0:2*L;
h = sin(2*pi*fm*(k-L))./(pi*(k-L));
h(L+1) = 2*fm;

%Change impulse response to frequency response
H = fft([h, zeros(1,200)]);
f = linspace(-1/2, 1/2, length(H));

%Plot frequency response
plot(f, [abs(fftshift(H))]); 
axis([-0.5 0.5 0 1.25]);
xlabel('$\hat f$', 'Interpreter', 'latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title('Truncated Sinc Frequency Response');

a_raw = x.*exp(-1i*2*pi*(1:N)'*(f_i - f_LO)/Fs);
a_raw_X = fftshift(fft(a_raw));
sinc_filtered = conv(real(a_raw), h);

%soundsc(real(sinc_filtered), Fs);
%sinc_filtered_2 = conv(real(sinc_filtered), h);
%soundsc(real(sinc_filtered_2), Fs);
%There is still a high pitched sound in the background, but when convoluted
%twice with the filter it gets better.

%% 4.4.2(b) Using FIR filters (Parks-Mclellan filter)

delta = 1/40;
h_pm = firpm(2*L, [0 fm fm+delta 0.5]*2, [1 1 0 0]);

%frequency response
H_pm = fft([h_pm, zeros(1,200)]);
f = linspace(-1/2, 1/2, length(H_pm));

%plot frequency response
plot(f,abs(fftshift(H_pm)));
xlabel('$\hat f$','Interpreter','latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title('Parks-Mclellan frequency response');

a_raw = x.*exp(-1i*2*pi*(1:N)'*(f_i - f_LO)/Fs);
a_raw_X = fftshift(fft(a_raw));
hmc_filtered = conv(real(a_raw), h_pm);

%soundsc(real(sinc_filtered), Fs);
%hmc_filtered_2 = conv(real(hmc_filtered), h_pm);
%soundsc(real(hmc_filtered_2), Fs);

%% Questions 4.2-4.6

%Question 4.2: The cutoff frequency is where the frequencies get destroyed
%by the lowpass filter. If this is set too low, some of the radio frequency
%may be cut off. If it is too high, the noise from other stations can e
%heard. The frequency response of the filter wil lhave a wider width if the
%cutoff frequency is larger, and vise versa. It can be shown by the
%frequency response where the singals are filtered out.

%Question 4.3: The filter length affects the quality by making the
%frequency response have more level magnitudes (shorter slide down). This
%can lower some of the noise due to the other stations.

%Question 4.4: The truncated sinc has better high frequency suppression.
%You can see this by the badwidth of the sinc being smaller than the
%Parks-Mclellan filter. Although the truncated sinc has a larger magnitude
%for the coefficients at the edges of the bandwideth, they are still at
%lower frequencies. 

%Question 4.5: The truncated sinc has better sound quality because the
%background high pitched noise is less evident than the Parks-Mclellan
%filter.

%Question 4.6: The casecade of the truncated sinc filters improved the
%quality by a good amount, but the Parks-Mclellan was less evident. The
%Parks-Mclellan model has magnitudes that are all basically 1 or 0, so
%multiplying these does not change much. The truncated sinc on the other
%hand has a wider range of magnitudes so this multiplication by cascading
%them actually improves teh sound quality.

%% 4.5 Listening to the Stations

clear;
load shortwave.mat
f_LO = 6e6;
x = raw(:,1) + 1i*raw(:,2);
X = fftshift(fft(x));
N = length(X);
bandwidth = 4000;
freqs = linspace(f_LO - Fs/2, f_LO + Fs/2, N);

%Change the range to find the different stations
minimum = 6.08e6;
maximum = 6.09e6;

%Find the indexes
minInd = min(find(freqs <= maximum & freqs >= minimum));
maxInd = max(find(freqs <= maximum & freqs >= minimum));
[maxVal, ind0] = max(abs(X(minInd:maxInd)));
f_i = freqs(minInd+ind0)
a_raw = x.*exp(-1i*2*pi*(1:N)'*(f_i - f_LO)/Fs);
a_raw_X = fftshift(fft(a_raw));
fhat = linspace(-1/2, 1/2, N);
new_X = a_raw_X.*(abs(fhat') <= bandwidth/Fs);

%soundsc(abs(ifft(new_X)), Fs);

%% Questions 4.7-4.8

%Question 4.7: The British man said there are 60,000 developers around
%China, found at 6.0201 MHz.

%Question 4.8: The low interest finance open house for HitchHouse Motor
%home show and sale, found at 6.0700 MHz.