% fMRI 
% 
clear
clf

% load 122 x 1 reference signal waveform 

load fmrisig.mat
s = s-mean(s); % remove average (DC) value
figure(1)
plot(s)
title('Reference signal')

% plot magnitude of fft of reference

figure(2)
plot(abs(fft(s)))
title('Magnitude of DFT of reference signal')

% load fMRI data 

load fmri.mat

% contains 122 x 64 x 64 array called 'x', consisting of 122 MRI 
% brain images of 64 by 64 pixels collected at 4 Hz  sampling rate

% each image is 64x64 pixels 
% time-length of image sequence = 122


% display MRI at time 1

figure(3)
im = reshape(x(1,:,:),64,64);
image(im/50)  % the divide by 50 scales the colormap to a nice range
title('Brain MRI at time t=1')

% plot the time series for pixel at (32,32) as a function of time

figure(4)
y = x(:,32,32);
plot(y)
title('Pixel (32,32) time course')

% plot the magnitude of the DFT for pixel at (32,32)

yz = y - mean(y); % remove DC value
figure(5)
plot(abs(fft(yz)))
title('Pixel (32,32) DFT magnitude')
