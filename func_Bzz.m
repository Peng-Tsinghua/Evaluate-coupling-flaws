function Bzz=func_Bzz(theta,V,B)
n=length(V);
A=zeros(n,n);
D=A;C=A; %[A,D';D,C]
for i=1:n
    C(i,i)=-B(i,i);
    for j=[1:i-1,i+1:n]
        C(i,j)=-B(i,j)*cos(theta(i)-theta(j)); % \p v_j \p v_i
        A(i,j)=-B(i,j)*V(i)*V(j)*cos(theta(i)-theta(j)); % \p theta_j \p theta_i
        D(j,i)=-B(i,j)*V(j)*sin(theta(i)-theta(j)); % \p theta_j \p v_i
        A(i,i)=A(i,i)-A(i,j);
        D(i,i)=D(i,i)+B(i,j)*V(j)*sin(theta(i)-theta(j));
    end
end
Bzz=[A,D;D',C];