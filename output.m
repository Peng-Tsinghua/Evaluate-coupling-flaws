clear;clc
Ml=zeros(10,1);SDl=Ml;MaxN=Ml;MinN=Ml;C1=Ml;C2=Ml;
%%
for Case=1:10
    %%
    switch Case
    case 1
        load('0IEEE9bus');
    case 2
        load('0IEEE14bus');
    case 3
        load('0IEEE39bus');
    case 4
        load('0IEEE118bus');
    case 5
        load('0IEEE145bus');
    case 6
        load('0Illinois200');
    case 7
        load('0IEEE300bus');
    case 8
        load('0Sg500');
    case 9
        load('0RTE1888');
    case 10
        load('0Polish3120sp');
    otherwise
        disp('Case wrong!')
    end
    %% output system dimension
%     % check phase cohesiveness
%     Bbus=(Bbus+Bbus')./2;%对称消除误差
%     G=graph(Bbus-diag(diag(Bbus)));%出掉self-loop
%     E=incidence(G);
%     disp([Case,min(min(cos(E'*angle(Vd))))])
%     %% output mean and var
%     disp([Case,mean(l),std(l),min(nn),max(nn)])
%     Ml(Case)=mean(l);SDl(Case)=std(l);MaxN(Case)=max(nn);
%     MinN(Case)=min(nn);
    %% check condition 1
    % 结果：所有cases都满足
    ct=zeros(1000,1);
    % Bbus0: whiteout bi
    Bbus0=Bbus-diag(diag(Bbus));
    Bbus0=Bbus0-diag(sum(Bbus0));
    for i=1:1000
        U=Vd(:,i);
        I=1i*Bbus0*U;
        S=U.*conj(I);
        P=real(S);Q=imag(S);
        V=abs(U);
        shift=0:0.01:pi;
        for j=1:length(shift)
            the=angle(U)+shift(j);
            c2t=sum(2*P.*cos(2*the)./V.^2-Q.*(1-sin(2*the))./V.^2);
            ct(i)=c2t;
            if ct(i)<0
                break
            end
        end
    end
    disp([Case,max(ct)]);
    C1(Case)=max(ct);
    %% Check condition 2
    % 结果：case5和case7不满足，其他都满足
    Bbus=(Bbus+Bbus')./2;%对称消除误差
    G=graph(Bbus-diag(diag(Bbus)));%出掉self-loop
    E=full(incidence(G));
    
    sum_the=zeros(1000,1);
    for i=1:1000
        U=Vd(:,i);
        V=abs(U);
        shift=0:0.01:pi;
        for j=1:length(shift)
            the=angle(U)+shift(j);
            temp=max(mod(abs(E')*the,2*pi));
            sum_the(i)=temp;
            if temp<pi
                break
            end
        end
    end
    disp([Case,max(sum_the)]);
    C2(Case)=max(sum_the);

%     %% Any possible correlation between \lambda and |theta-psi|?
%     % result: not very good, even after the degree constrain
%     maxthpsi=zeros(1000,1);
%     Bbus=(Bbus+Bbus')./2;%对称消除误差
%     G=graph(Bbus-diag(diag(Bbus)));%出掉self-loop
%     deg=degree(G);
%     mdeg=mean(deg);
%     index=find(deg>mdeg);
%     n=length(index);
%     for i=1:1000
%         U=Vd(:,i);
%         the=angle(U);
%         psi=zeros(n,1);
%         for j=1:n
%             psi(j)=angle(sum(U([1:index(j),index(j)+1:end]))/(n-1));
%         end
%         maxthpsi(i)=max(abs(the(index)-psi));
%     end
%     figure, scatter(maxthpsi,l)
end