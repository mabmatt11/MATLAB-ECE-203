function [xx, tt] = one_cos(amp, ff, phase, dur)
tt = 0:(2*pi)/(20*ff):dur;
xx = amp*cos(ff*tt + phase);
