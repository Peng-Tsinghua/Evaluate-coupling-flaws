function [E,del]=Vthe2Edel(V,the,xdd,P,Q)
% 将机端电压V,the转换为发电机内电势E,del
% P,Q为机端输出功率;xdd: xd'
I=conj((P+1i*Q)./(V.*exp(1i*the)));
U=V.*exp(1i*the)+I.*1i.*xdd;
E=abs(U);
del=angle(U);