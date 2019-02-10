function y = ADSR(t,t1,t2,t3,t4,A1,A2,A3)
% Sum of 4 straight lines:
% line 1: (0,0) to (t1,A1)
% line 2: (t1,A1) to (t2,A2)
% line 3: (t2,A2) to (t3,A3)
% line 4: (t3,A3) to (t4,0)
y = myline(t,0,0,t1,A1) + myline(t,t1,A1,t2,A2) + ...
myline(t,t2,A2,t3,A3) + myline(t,t3,A3,t4,0);