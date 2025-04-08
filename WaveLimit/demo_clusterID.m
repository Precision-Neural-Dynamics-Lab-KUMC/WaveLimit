
nd = 20;
d = -5:0.1:5;
sigma = 1;
sigma2 = 2.5;

x = sort(3*randn(1000,nd));

y = exp(-0.5 * ((  sum(x,2)  )./sigma).^2) ./ (sqrt( (2*pi)^nd .* sigma.^2));
y2 = exp(-0.5 * ((  sum(x,2)  )./sigma2).^2) ./ (sqrt( (2*pi)^nd .* sigma2.^2));
y(y2==0) = 0;
y2(y2==0) = eps;
sum(y./y2)

figure; plot(x,y)


d = -5:0.1:5
y = exp(-0.5*(x./1).^2);
figure; plot(x,y)
y2 = exp(-0.5*(x./2.5).^2);
hold on; plot(x,y2)