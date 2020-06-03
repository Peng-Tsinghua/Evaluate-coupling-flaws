function [E,del]=Vthe2Edel(V,the,xdd,P,Q)
% convert the terminal voltage phasor V,the to the inner voltage phasor E,del behind the transient reactance
% P,Q is the terminal output active and reactive power; xdd is the transient reactance xd'
I=conj((P+1i*Q)./(V.*exp(1i*the)));
U=V.*exp(1i*the)+I.*1i.*xdd;
E=abs(U);
del=angle(U);