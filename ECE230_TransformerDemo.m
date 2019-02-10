% ECE 230 exercise 15A transformer review
k = 1
XL1 = 100
n = 2;
R = 100;
Vs = 10;

%ideal transformer with R2
i2_ideal = Vs * n / R

%linear transformer with R2
XL2=n^2*XL1;    %inductance proportional to N^2
XM=k*sqrt(XL1*XL2);
A=[i*XL1 -i*XM;-i*XM (R + i*XL2)];
B=[Vs ; 0];
I=inv(A)*B;
i2_linear = abs(I(2))

%linear transformer without R2
i1_open_circuit_L2 = Vs/XL1