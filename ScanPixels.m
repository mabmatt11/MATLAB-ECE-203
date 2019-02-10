function Decision = ScanPixels(fmriData, k, MagnitudeFraction, MaxFreqIndexDiff)
% k(1) and k(2) are the indexes of the peaks we are looking for
maxval = max(fmriData(:));
[~,xmax,ymax] = size(fmriData);
Decision = false(xmax,ymax);
for i=1:xmax
   for j=1:ymax
      TimeSeries = fmriData(:,i,j);
      MagSpectrum = abs(fft(TimeSeries));
      MagSpectrum(1) = 0;                 % Remove DC
      maxMagSpectrum = max(MagSpectrum);
      l = find( MagSpectrum == maxMagSpectrum);
      Decision(i,j) = maxMagSpectrum >= MagnitudeFraction*maxval && abs(l(1)-k(1))<=MaxFreqIndexDiff && abs(l(2)-k(2))<=MaxFreqIndexDiff;
   end
end


