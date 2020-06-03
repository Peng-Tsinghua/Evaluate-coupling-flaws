% This function extends the original Y matrix to include xd'=xdd
% index(column): The indices of the nodes in the original matrix that should be extent. Each node will have an extra branch and node.
% y(column): impedence to be extent, p.u.
% The index in the extent node is next to the last node of the original system
function Ynew=Yextend(Yold,index,y)
n=size(Yold,1); % size of original network
m=length(y); % number of SG
index2=(1:m)+n;
Ytemp=sparse(n+m,n+m);
Ytemp(1:n,1:n)=Yold;
Ynew=Ytemp+sparse(index,index2,-y,n+m,n+m)+sparse(index2,index,-y,n+m,n+m)+sparse(index,index,y,n+m,n+m)...
    +sparse(index2,index2,y,n+m,n+m);
