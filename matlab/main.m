clear;clc;clf;

n = 1000; % Length of array (>30 is recommended).
r = 1.0; % Roughness (1 for omega-squared model).
x = StochasticSTF(n);
plot(x)

function y = StochasticSTF(n,r,d)
    if (~exist('r','var'))
         r = 1.0;
    end
    if (~exist('d','var'))
         d = 2;
    end
    n1 = floor(n*r/(1+r));
    n2 = n-n1+1;
    x1 = BesselBridge(n1+1,d);
    x2 = BesselBridge(n2+1,d);
     y = conv(x1(2:end),x2(1:end-1));
     y = y/sum(y);
end

function y = BesselBridge(n,d)
    x = BrownBridges(n,d);
    y = rms(x,2);
end

function y = BrownBridges(n,d)
    x = [zeros(1,d); cumsum(randn(n-1,d),1)];
    y = x-(0:n-1)'*x(n,:)/(n-1);
end
