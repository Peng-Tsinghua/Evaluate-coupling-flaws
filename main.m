% reference
% [1]. P. M. Anderson and A. A. Fouad, Power system control and stability. IEEE, 2003.
% [2] M. A. Pai, Energy Function Analysis for Power System Stability,Kluwer Academic Publishers, Boston, 1989.
% [3] Athay, T.; Podmore, R.; Virmani, S., "A Practical Method for the Direct Analysis of Transient Stability," IEEE Transactions on Power Apparatus and Systems , vol.PAS-98, no.2, pp.573-584, March 1979.
% [4] V. Vittal, “Transient stability test systems for direct stability methods,” IEEE Transactions on Power Systems, vol. 7, no. 1, pp. 37–43, Feb. 1992, doi: 10.1109/59.141684.
%% Load case
clear;clc;
define_constants;
mpopt=mpoption('pf.nr.max_it',50);

for Case=4 % 1--9bus; 2--39bus
switch Case
    case 1
        mpc0=loadcase('case9');
%         load('0IEEE9bus');
    case 2
        mpc0=loadcase('case14');
%         load('0IEEE39bus');
    case 3
        mpc0=loadcase('case39');
%         load('0IEEE118bus');
    case 4
        mpc0=loadcase('case118');
%         load('0IEEE145bus');
    case 5
        mpc0=loadcase('case145');
%         load('0Polish3120sp');
    case 6
        mpc0=loadcase('case_illinois200');
%         load('0IEEE14bus');
    case 7
        mpc0=loadcase('case300');
        mpc0=ext2int(mpc0);
%         load('0IEEE300bus');
    case 8
        mpc0=loadcase('case_ACTIVSg500');
        mpc0=ext2int(mpc0);
%         load('0Illinois200');
    case 9
        mpc0=loadcase('case1888rte');
        mpc0=ext2int(mpc0);
%         load('0RTE1888');
    case 10
        mpc0=loadcase('case3120sp');
%         load('0Sg500');
    otherwise
        disp('Case wrong!')
end
LDbus=find(mpc0.bus(:,PD)>0); % find load bus index
nLD=length(LDbus); % total number of loads
indexSG=find(mpc0.gen(:,GEN_STATUS)==1); % index of turn-on SG
indexGen=find(mpc0.gen(:,PG)>0); % index of SG in matrix mpc.gen in MATPOWER case. The rest are condensors
% SGbus=mpc0.gen(indexGen,1); % bus index of generators
SGbus=mpc0.gen(indexSG,1); 
nSG=length(SGbus); % total number of generators
n0=size(mpc0.bus,1); % total number of buses

switch Case
    case 1
        % 3 generators' parameters from ref. [1]
        % all reactance values are in per unit on a base 100MVA
        % all time constants are in s.
        xd=[0.146;0.8958;1.3125];
        xdd=[0.0608;0.1198;0.1813]; % xd'
        Td0=[8.96;6;5.89]; % T'd0
        H=[23.64;6.40;3.01]; % the inertia constants converted to system Sb=100MVA
    case 3
        % 10 generators' parameters from ref. [2]
        % all reactance values are in per unit on a base 100MVA
        % all time constants are in s. 
        % parameters are given in bus order 30-39, which corresponds to SG order 10,2-9,1
        % in ref. [2]
        xd=[.1;.295;.2495;.262;.67;.254;.295;.290;.2106;.02];
        xdd=[.031;.0697;.0531;.0436;.132;.05;.049;.057;.057;.006]; % xd'
        Td0=[10.2;6.56;5.7;5.69;5.4;7.3;5.66;6.7;4.79;7]; % T'd0
        H=[42;30.3;35.8;28.6;26.0;34.8;26.4;24.3;34.5;500]; % the inertia constants converted to system Sb=100MVA
    case 4
        % 20 generators' parameters from ref. [3]
        % rest condensors' xdd are estimated
        % all reactance values are in per unit on a base 100MVA
        S0=sqrt(mpc0.gen(:,PMAX).^2+mpc0.gen(:,QMAX).^2);
        Q0=mpc0.gen(:,QMAX);
        V0=mpc0.gen(:,VM);
        QV0=abs(Q0)./V0;
        xdd=zeros(nSG,1);    
        xdd(indexGen)=[0.0700    0.1167    0.0875    0.0438    0.1000    0.1400...
            0.0875    0.0875    0.1167    0.0875    0.0875    0.0700...
            0.0500    0.0636    0.0467    0.0538    0.0875    0.1750    0.0467]';
        % in this case, no off-state SG
        indexCon=1:nSG;
        indexCon(indexGen)=[];
        xdd(indexCon)=17.91./QV0(indexCon);
    case 5
        % 50 generators' parameters from ref. [4]
        % all reactance values are in per unit on a base 100MVA
        % parameters are given in bus order
        xdd=[0.4769;0.0213;0.1292;0.6648;0.5291;0.0585;1.6000;...
            0.3718;0.0240;0.0839;0.1619;0.4824;0.2125;0.0795;0.1146;...
            0.1386;0.0924;0.0135;0.1063;0.0122;0.0208;0.0312;0.0248;...
            0.2029;0.0240;0.0122;0.0924;0.0024;0.0022;0.0017;0.0014;...
            0.0002;0.0017;0.0089;0.0017;0.0001;0.0010;0.0001;0.0016;...
            0.0003;0.0008;0.0001;0.0004;0.0001;0.0003;0.0001;0.0003;...
            0.0023;0.0004;0.0018];

    otherwise
        % get nominal operation point by AC-pf solover
        result0=runpf(mpc0);    
        % when Pmax or Qmax is not provided, we use the nominal output to
        % replace it
        S0=sqrt(result0.gen(indexSG,QG).^2+result0.gen(indexSG,PG).^2);
        Q0=max(mpc0.gen(indexSG,QMAX),1.2*result0.gen(indexSG,QG));
        V0=mpc0.gen(indexSG,VM);
        QV0=abs(Q0)./V0;
        xdd=17.91./QV0;
end
%% solve nominal case
% for each generator xd' is added to the network
% all passive nodes are reduced via Kron reduction
Yold=makeYbus(mpc0);
Ypre=Yextend(Yold,SGbus,-1i./xdd);
index_left=[LDbus;(n0+1:n0+nSG)']; % order: load bus + SG bus
n=length(index_left); % bus number of reduced network
index_reduce=1:n0+nSG;
index_reduce(index_left)=[];
Yre=KronRe_new(Ypre,index_reduce',index_left);
Bbus=full(imag(Yre));
Bbus=1/2*(Bbus'+Bbus); % 对称化消除计算误差
% get nominal operation point by AC-pf solover
result0=runpf(mpc0);
V0=result0.bus(:,VM);the0=result0.bus(:,VA)/180*pi;
[E,del]=Vthe2Edel(V0(SGbus),the0(SGbus),xdd,...
    result0.gen(indexSG,PG)/mpc0.baseMVA,...
    result0.gen(indexSG,QG)/mpc0.baseMVA); %steady-state发电机内电势
V=[V0(LDbus);E]; the=[the0(LDbus);del];
U=V.*exp(1i*the);
I=1i*Bbus*U;
S=U.*conj(I);
P=real(S);Q=imag(S);
% calculate \lambda in the nominal case
% Bbus0: whiteout bi
Bbus0=Bbus-diag(diag(Bbus));
G=graph(Bbus0);Incid=incidence(G);
Bbus0=Bbus0-diag(sum(Bbus0));


Bzz=func_Bzz(the,V,Bbus0);
Schur=Bzz(n+1:end,n+1:end)-Bzz(1:n,n+1:end)'*pinv(Bzz(1:n,1:n))*Bzz(1:n,n+1:end);
l0=min(eig((Schur+Schur')/2)); % 对称化消除计算误差
% lB=min(eig(-Bbus));
c1=-2*sum(sum(Bbus));% -2*sum(bi)
c20=sum(2*P.*cos(2*the)./V.^2-Q.*(1-sin(2*the))./V.^2);
c0=c1+c20;
% disp([Case,l0,lB,c1,c0])

N=1000; % number of all cases
l=zeros(N,1);c=l;c2=c;nn=l;
Vd=zeros(n,N);
% the nominal case as the 1st case
Vd(:,1)=U;
l(1)=l0;
nn(1)=length(find(eig(Schur)<0)); 
c2(1)=c20;
c(1)=c0;

%% random scenario
i=1; % initialize the counter
sigma=0.3; % standard deviation for Gaussian distributions
while i<N
    % random volatile smart grid scenario
    mpc1=mpc0;
    % fluactuating loads 50%
    mpc1.bus(LDbus,PD)=mpc0.bus(LDbus,PD)+0.5*sigma*randn(nLD,1).*mpc0.bus(LDbus,PD);
    
    % renewables with stochastic power generation 30%
    mpc1.gen(indexGen,PG)=mpc0.gen(indexGen,PG)+0.3*sigma*randn(length(indexGen),1).*mpc0.gen(indexGen,PG);
    
    % dispatch imbalanced power to fast-ramping gneration 20% and controllable loads 5%
    Pim=sum(mpc1.bus(LDbus,PD)-mpc0.bus(LDbus,PD))-sum(mpc1.gen(indexGen,PG)-mpc0.gen(indexGen,PG)); % imbalanced power: increased load - increased generation
    Pca=0.2*sum(result0.gen(indexGen,PG))+0.05*sum(mpc0.bus(LDbus,PD)); % dispathable capacity
    mpc1.gen(indexGen,PG)=mpc1.gen(indexGen,PG)+Pim*0.2.*result0.gen(indexGen,PG)./Pca; % dispatch fast-ramping generation
    mpc1.bus(LDbus,PD)=mpc1.bus(LDbus,PD)-Pim*0.05.*mpc0.bus(LDbus,PD)./Pca; % dispatch controllable loads
    mpc1.bus(LDbus,QD)=mpc0.bus(LDbus,QD).*mpc1.bus(LDbus,PD)./mpc0.bus(LDbus,PD); % assume constant power factor
    % calculate \lambda
    % try to solve volatile case by AC-pf solover
    result=runpf(mpc1,mpopt);
    % if pf converged with cohesive phase, then get the operation point and do calculations
    if result.success
        V0=result.bus(:,VM);the0=result.bus(:,VA)/180*pi;
        [E,del]=Vthe2Edel(V0(SGbus),the0(SGbus),xdd,...
            result0.gen(indexSG,PG)/mpc0.baseMVA,...
            result0.gen(indexSG,QG)/mpc0.baseMVA); %steady-state发电机内电势
        V=[V0(LDbus);E]; the=[the0(LDbus);del];
        if min(cos(Incid'*the))>=0 % check cohesive phase
            i=i+1;
            U=V.*exp(1i*the);
            I=1i*Bbus*U;
            S=U.*conj(I);
            P=real(S);Q=imag(S);
            % calculate \lambda in the nominal case
            Bzz=func_Bzz(the,V,Bbus0);
            Schur=Bzz(n+1:end,n+1:end)-Bzz(1:n,n+1:end)'*pinv(Bzz(1:n,1:n))*Bzz(1:n,n+1:end);
            nn(i)=length(find(eig((Schur+Schur')/2)<0)); % number of negative eigenvalues
            l(i)=min(eig((Schur+Schur')/2));
            c2(i)=sum(2*P.*cos(2*the)./V.^2-Q.*(1-sin(2*the))./V.^2);
            c(i)=c1+c2(i);
            Vd(:,i)=U;
        end
    end
end
%%
figure,plot(l,'LineStyle','none','Marker','.');
ylabel('\lambda')
% figure,plot(c,'LineStyle','none','Marker','.');
% ylabel('C')
%% Save data
switch Case
    case 1
        save('0IEEE9bus','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 2
        save('0IEEE14bus','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 3
        save('0IEEE39bus','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 4
        save('0IEEE118bus','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 5
        save('0IEEE145bus','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 6
        save('0Illinois200','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 7
        save('0IEEE300bus','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 8
        save('0Sg500','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 9
        save('0RTE1888','Bbus','l0','c1','Vd','l','nn','c2','c');
    case 10
        save('0Polish3120sp','Bbus','l0','c1','Vd','l','nn','c2','c');
end
end