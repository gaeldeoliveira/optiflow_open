t = (1:0.1:5)';

% model
px = [ 5 2 1 ];
x =  polyval(px,t);

py = [ -2 1 1 ];
y = polyval(py,t);

pz = [ 1 20 1 ];
z = polyval(pz,t);

%  plot model
figure
plot3(x,y,z)
hold all

% simulate measurement 
xMeasured = x+2*(rand(length(x),1)-0.5);
yMeasured = y+2*(rand(length(y),1)-0.5);
zMeasured = z+2*(rand(length(z),1)-0.5);

% plot simulated measurements
plot3(xMeasured, yMeasured, zMeasured,'or')
hold off
grid on

% least squares fit 
A = [t.^2, t, t./t];
pxEstimated = A\xMeasured;
pyEstimated = A\yMeasured;
pzEstimated = A\zMeasured;