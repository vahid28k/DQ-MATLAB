function [c,cdm]=quad_int_mul_u_sens(a,b)
d=size(b,2);cmatrix=[];cdmatrix=[];
for i=1:d
    %%%change it to pol_mul_jacobi for Legendre and Chebyshev
    [ci, cdi]=pol_mul_g(a(i),b(:,i));
    cmatrix=[cmatrix ci];
    cdmatrix=[cdmatrix cdi];
end;
c=ones(size(b,1),1);
for i=1:d
    c=c.*cmatrix(:,i);
end;
cdm=zeros(size(b,1),d);
for r=1:d
ind=[1:d];ind(r)=[];
cd=ones(size(b,1),1);
for i=1:(d-1)
    cd=cd.*cmatrix(:,ind(i));
end;
cd=cd.*cdmatrix(:,r);
cdm(:,r)=cd;
end;



%rr=[c1 c2 c3 c4 c5].... c*(cd1/c1)