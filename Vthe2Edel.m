function [E,del]=Vthe2Edel(V,the,xdd,P,Q)
% �����˵�ѹV,theת��Ϊ������ڵ���E,del
% P,QΪ�����������;xdd: xd'
I=conj((P+1i*Q)./(V.*exp(1i*the)));
U=V.*exp(1i*the)+I.*1i.*xdd;
E=abs(U);
del=angle(U);