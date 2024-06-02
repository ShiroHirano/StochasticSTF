clear;clc;clf;

n = 1000; % Length of array (>30 is recommended).
r = 1.0; % Roughness (1 for omega-squared model).
x = StochasticSTF(n,r);
plot(x)

function out=StochasticSTF(n,r,d)
    if (~exist('d','var'))
         d = 2
    end
    n1=floor(n*r/(1+r));
    n2=n-n1+1;
    x1=BesselBridge(n1+1,d);
    x2=BesselBridge(n2+1,d);
    out=conv(x1(2:end),x2(1:end-1));
    out=out/sum(out);
end

function out=BesselBridge(n,d)
    x=BrownBridges(n,d);
    out=rms(x,2);
end

function out=BrownBridges(n,d)
    x=[zeros(1,d); cumsum(randn(n-1,d),1)];
    out=x - (0:n-1)'*x(n,:)/(n-1);
end
