function [c cd]=cons_computation(x,par)
%%Instead of enforcing constraint -1<g<1 we enforce more conservative
%%constraints i.e. -1+delta<g<1-delta 
delta=1e-6;

%%for Gaussian make the limit high
lim=1e5-delta;
c=heaviside(abs(x)-lim).*(abs(x)-lim).^2;c=c*par;
cd=heaviside(abs(x)-lim).*(2*(abs(x)-lim)).*sign(x);cd=cd*par;

