function y = envelope(tin,seconds_per_pulse)
% Shift time from [tin(1),tin(end)] to [0,T],
% where T = tin(end)-tin(1).
t = tin - tin(1);
T = t(end); % Remember that we have already just subtracted
% tin(1) from all the elements of tin.
% No matter how long the note is, the attack, the delay, and
% the release are always of fixed duration. Instead of
% fixing the duration of the sustain, we fix its slope so
% that the sustain would connect the end of the delay to the
% beginning of the release when the delay ends at (t2,A2)
% and the release starts at (t3,A3) in the case that
% t2 = .2 and t3 = .8. These choices for t2 and t3 were
% based on plots of the envelope when [0,T] = [0,1].
% It remains to compute the value of the sustain at the
% true value of t3. We also prevent the sustain from
% becoming negative, which could otherwise happen in a note
% of very long duration.
t1 = 1/8 * seconds_per_pulse;
A1 = 1;
t2 = .2 * seconds_per_pulse;
A2 = .85;
t3 = T - .2 * seconds_per_pulse;
A3 = .75;
t4 = T;
% When t2 = .2 and t3 = .8, t3-t2 = .6,
% giving us the slope of the sustain.
A3p = (A3-A2)/.6*(t3-t2) + A2;
if A3p < 0 % prevent sustain from becoming negative
    t3 = t2 + .6*A2/(A2-A3);
    A3p = 0;
end
y = ADSR(t,t1,t2,t3,t4,A1,A2,A3p);