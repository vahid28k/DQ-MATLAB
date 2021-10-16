function [c cd]=pol_mul_g(a,b)

if a==0
    c = ones(length(b),1);
    cd = zeros(length(b),1);
end;
if a==1
    c = b;
    cd = ones(length(b),1);
end;

if a>1
i=2;
P(:,1) = ones(length(b),1);
Pd(:,1) = zeros(length(b),1);
P(:,2) = b;
Pd(:,2) = ones(length(b),1);
while i<=a

P(:,i+1) = b.*P(:,i) - ((i-1)*P(:,i-1));
Pd(:,i+1) = (P(:,i)+b.*Pd(:,i)) - ((i-1))*Pd(:,i-1);
i=i+1;

end;
c = P(:,i);
cd = Pd(:,i);
end;


