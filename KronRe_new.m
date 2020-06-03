% Kron reduction 
% re: nodes indices to be reduced, column
% le: nodes indices to be left,ordered, column
function Ynew=KronRe_new(Yold,re,le)
% reorder Yold
order=[le;re];
Yold=Yold(order,order);
m=size(Yold,1); % size of original network
% reduction: left the first n nodes
n=length(le);
Yrr=Yold(n+1:m,n+1:m);
Yrn=Yold(n+1:m,1:n);
Ynr=Yold(1:n,n+1:m); % dont use transpose!! for image it's conj transpose!
Ynn=Yold(1:n,1:n);
Ynew=Ynn-Ynr/Yrr*Yrn;