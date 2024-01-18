esno = 10^((5)/10);

x =0:0.1:2;
y = 1./(x.^2);
z = exp(-esno*((x.^2)/4));
figure
plot(x,y)
hold on;
plot(x,z)
grid