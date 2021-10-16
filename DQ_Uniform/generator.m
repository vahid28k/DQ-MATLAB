function [XW,deltamain]=generator(d,p,n_s)

%%number of terms for total order
n_terms=floor(factorial(d+p)/(factorial(d)*factorial(p)));
aind = total_degree_indices(d, p);

%%%initialization
b=lhsdesign(n_s,d)*2-1;
w=(n_terms/n_s)*ones(n_s,1);

delta=1;count=1;
xnew=[b(:,1)];
for i=2:d
    xnew=[xnew;b(:,i)];
    
end;
xnew=[xnew;w];

RHS=zeros(n_terms+n_s*(d+1),1);RHS(1,1)=1;


%%%%The tolerance is set to \epsilon=1e-9.
while delta>10^(-9)
    
    xold=xnew;
    
    for i=1:d
        b(1:n_s,i)=xold(1+(i-1)*n_s:i*n_s,1);
    end;
    w=xold(d*n_s+1:(d+1)*n_s,1);
    
    
    R=zeros(n_terms,1);J=zeros(n_terms,(d+1)*n_s);
    for i=1:n_terms
        a=aind(i,:);
        %%The function pol_mul_jacobi inside quad_int_mul_u_sens computes the
        %%jacobi polynomials. For uniform weight alpha=0, beta=0 is pre-set.
        [Ri,Rdij]=quad_int_mul_u_sens(a,b);
        R(i,1)=w'*Ri;
        for j=1:d
            Rsens=Rdij(:,j);
            J(i,1+(j-1)*n_s:j*n_s)=w'.*Rsens';
        end;
        J(i,d*n_s+1:(d+1)*n_s)=Ri';
    end;
    
    
    %%%computing residual and its jacobian for constraints
    Rcons=[];Jcons=[];
    for i=1:(d)
        parm = 1/norm(R-[n_terms;zeros(n_terms-1,1)]);
        parm=max(parm,1000);
        [res,resd]=cons_computation(b(:,i),parm);
        Rcons(1+(i-1)*n_s:i*n_s,1)= res;
        Jcons(1+(i-1)*n_s:i*n_s,1)= resd;
    end;
    [resw,resdw]=cons_w_computation(w,parm);
    Rcons=[Rcons;resw];
    Jcons=[Jcons;resdw];
    
    %%Putting together the augemented residual and its jacobian
    Jcons=[diag(Jcons)];
    R=[R;Rcons];R=R-RHS;
    J=[J;Jcons];
    
    
    
    delta=norm(R);
    [z1,z2,z3]=svd(J);
    
    %%%Selecting regularization parameter based on the residual norm. These are default values that can be changed
    %%%in the case of e.g. higher order p and dimension d for better convergence.
    if delta<2000 dtikh=1000; end;
    if delta<500 dtikh=500; end;
    if delta<200 dtikh=100; end;
    if delta<50 dtikh=30; end;
    if delta<10 dtikh=10; end;
    if delta<1 dtikh=5; end;
    if delta<0.5 dtikh=1; end;
    if delta<0.1 dtikh=0.1; end;
    if delta<0.01 dtikh=0.05; end;
    if delta<0.001 dtikh=0.01; end;
    if delta<0.0001 dtikh=0.001; end;
    if delta<10^(-5) dtikh=5*10^(-4); end;
    
    %%%Filtering the small singular values.
    z2inv=zeros((d+1)*n_s,n_terms+(d+1)*n_s);
    for i=1:min((d+1)*n_s,n_terms+(d+1)*n_s)
        %  if (  z2(i,i)>0 &&   delta>1e-3 )
        %      ff=(z2(i,i)^2)/(z2(i,i)^2+dtikh^2);
        %      z2inv(i,i)=ff/(z2(i,i));
        %  else
        z2inv(i,i)=1/(z2(i,i)+dtikh);
        %  end;
    end;
    
    %%computing Newton's step
    def=(z3*z2inv*z1')*R;
    xnew = xold - def;
    
    %%computing Newton decrement to assess convergence
    ndecr=sqrt(def'*(J'*R));
    
    %%%Saving some variables and showing the residual norm
    count=count+1;
    if mod(count,1)==0
        res=delta
        iters(count)=delta;
        tikhcons(count)=dtikh;
        ndef(count)=norm(def);nobj(count)=norm(J*def-R);ndcr(count)=ndecr;
        sol(:,:,count)=b;
    end;
    
    
    deltamain=delta;
    
    
    if ((delta/ndecr)>10000)
        msg = ['[Newton Decrement / ||R||] is so small: ',num2str(ndecr/delta),'. Increase the number of points.'];
        disp(msg);
        break;
    end;
    
    ccount=5000;
    if (count>ccount)
        msg = ['Too many iterations: ',num2str(ccount)];
        disp(msg);
        break;
    end;
    
    if delta>10^30
        msg = ['Failed miserably! The regularization parameter is too small'];
        disp(msg);
        break;
    end;
    
    
end;

XW=[b/2+0.5 w];


subplot(1,2,1)
scatter(xnew(1:n_s,1),xnew(n_s+1:2*n_s,1),'filled');axis square;xlim([-1 1]);ylim([-1 1]);grid on;xlabel('x^{(1)}');ylabel('x^{(2)}');
subplot(1,2,2);
semilogy(iters);xlabel('iteration');ylabel('$\| R \|$','interpreter','latex');axis square;grid on;
