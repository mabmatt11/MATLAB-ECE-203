h1 = [1/2, 1/2];
h2 = [1/2, 1/2];

for R = 1:96
    h1 = conv(h1, h2);
end
