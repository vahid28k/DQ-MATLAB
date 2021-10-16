function [c cd]=pol_mul_jacobi(a,b)
%%change alpha=-0.5;beta=-0.5 for Chebyshev polynomial.
alpha=0;beta=0;
if a==0
    c = ones(length(b),1);
    cd = zeros(length(b),1);
end;
if a==1
    c=[gamma(alpha+2)/gamma(alpha+beta+2)]*[(gamma(alpha+beta+2)/gamma(alpha+1)) + (gamma(alpha+beta+3)/gamma(alpha+2)) * ((b-1)/2)  ];  
    cd =0.5*ones(length(b),1)*[gamma(alpha+beta+3)/gamma(alpha+beta+2)];
end;

if a>1
i=2;
P(:,1) = ones(length(b),1);
Pd(:,1) = zeros(length(b),1);
P(:,2) = [gamma(alpha+2)/gamma(alpha+beta+2)]*[(gamma(alpha+beta+2)/gamma(alpha+1)) + (gamma(alpha+beta+3)/gamma(alpha+2)) * ((b-1)/2)  ]; 
Pd(:,2) = 0.5*ones(length(b),1)*[gamma(alpha+beta+3)/gamma(alpha+beta+2)];
while i<=a
    
cons_L1=2*i*(i+alpha+beta)*(2*i+alpha+beta-2);
cons_R1=(2*i+alpha+beta-1)*[(2*i+alpha+beta)*(2*i+alpha+beta-2)*b+alpha^2-beta^2];
cons_R2=2*(i+alpha-1)*(i+beta-1)*(2*i+alpha+beta);
P(:,i+1) = [cons_R1.*P(:,i) - (cons_R2*P(:,i-1))]/cons_L1;
Pd(:,i+1) = ([(2*i+alpha+beta-1)*(2*i+alpha+beta)*(2*i+alpha+beta-2)]*P(:,i)+(cons_R1).*Pd(:,i))/cons_L1 - (cons_R2/cons_L1)*Pd(:,i-1);


%hermite
% P(:,i+1) = b.*P(:,i) - ((i-1)*P(:,i-1));
% Pd(:,i+1) = (P(:,i)+b.*Pd(:,i)) - ((i-1))*Pd(:,i-1);


i=i+1;

end;
c = P(:,i);
cd = Pd(:,i);
end;