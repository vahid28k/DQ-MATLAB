clc;clear all;close all;

%%%Some successful designed quadratre rules that are already found reported below. Note that they confirm Stroud rule: n_s=d+1 for p=2 and n_s=2d for p=3.
%%d=2,p=[1:1:10];n_s=[1,3,4,6,7,10,12,16,17,23];
%%d=3,p=[1:1:10];n_s=[1,4,6,10,13,22,26,43,51,74];

%%%Some of these rules are found by testing many intial guesses i.e.
%%%several re-runs. If the same number of points is not found simply
%%%increase the number of points.

%%%First d columns of XW are the nodes and last column is the weights. 

d=2;p=10;n_s=23;

tic
[XW,deltamain]=generator(d,p,n_s);
toc


%%%For Gaussian e.x. check \int_{\R^25} x_1^2 + x_2^2 \rho(x_1,...,x_25) dxdy = 2
int_2_var=[XW(:,1).^2+XW(:,2).^2]'*XW(:,d+1);

%%%or square integral of all d dimensions = 25
int_25_var=sum(XW(:,1:d).^2,2)'*XW(:,d+1);