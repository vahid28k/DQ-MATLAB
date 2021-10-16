function [c cd]=cons_w_computation(x,par)
%%Instead of enforcing constraint 0<w we enforce more conservative
%%constraints i.e. delta<w  
delta=1e-6;
x=x-delta;
c=0.25*(x-abs(x)).^2;c=c*par;
cd=(0.5*(x-abs(x))).*(1-sign(x));cd=cd*par;

