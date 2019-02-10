function y = myline(t,t1,y1,t2,y2)
y = ( (y2-y1)/(t2-t1)*(t-t1) + y1 ) .* (t1<=t & t<t2);