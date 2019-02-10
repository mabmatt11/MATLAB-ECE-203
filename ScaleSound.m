function y = ScaleSound(x)
% Do the scaling operation performed by soundsc and
% return the result. Let a = min(x(:)), b = max(x(:)), and
%
% f(x) = 2(x-a)/(b-a) - 1.
%
% Notice that f is linear and satisfies
%
% f(a) = -1 and f(b) = +1.
a = min(x(:));
b = max(x(:));
if b==a % then x is a constant vector, and we put
    y = zeros(size(x));
else
    y = (x-a)/(b-a)*2-1;
end