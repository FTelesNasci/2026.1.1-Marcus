close all
clear all
clc

%% Article - johansson 2000 => The Quadruple Tank Process

% Example 2
A=[ -0.764  0.387 -12.9 0 0.952  6.05;
     0.024 -0.174  4.31 0 -1.76 -0.416;
     0.006 -0.999 -0.0578 0.0369 0.0092 -0.0012;
     1 0 0 0 0 0; 0 0 0 0 -10 0; 0 0 0 0 0 -5];
B=[0 0 0 0 20 0;0 0 0 0 0 10]'; 
C=[0 0 1 0 0 0;0 0 0 1 0 0]; D=zeros(2);

delta=0.6; % original
% delta=0.85 % 0.9995
Delta1=diag([1+delta 1+delta]); 
Delta2=diag([1-delta 1-delta]); 

% Modelo aumentado
% Sistema nominal
Ahn=[A zeros(size(B)); -C zeros(size(D))];
Bhn=[B; zeros(size(D))];
Chn=[C zeros(size(D))]; Dhn=D;

% Sistema incerto
Ah{1}=[A zeros(size(B)); -C zeros(size(D))];
Ah{2}=[A zeros(size(B)); -C zeros(size(D))];

Bh{1}=Bhn*Delta1; Bh{2}=Bhn*Delta2;

Ch{1}=[C zeros(size(D))];
Ch{2}=[C zeros(size(D))];
Dh{1}=D; Dh{2}=D;

%% Parâmetros do controlador
gamma=1;

G1=[eye(2); eye(2); eye(2); eye(2)];
G2=[eye(2); eye(2); eye(2); eye(2)];

Qc=diag([1 1 1 1 1 1 1e5 1e5]); Rc=eye(2);

[n m]=size(Bh{1});  % Tamanho da planta 

X1=sdpvar(n, n, 'symmetric');
X2=sdpvar(n, n, 'symmetric');

L1=sdpvar(m, n, 'full');
L2=sdpvar(m, n, 'full');

W1=sdpvar(m, m, 'symmetric');
W2=sdpvar(m, m, 'symmetric');


%% LMI's
LMIs=[ ];
LMIs=[LMIs, X1>=0, X2>=0];
LMIs=[LMIs, W1>=0, W2>=0];

%% Subsistema 01

Euclideana_1=(Ah{1}*X1+Bh{1}*L1)'+(Ah{1}*X1+Bh{1}*L1);

S_01=[ (Euclideana_1-gamma*X1) (gamma*X1)'       X1'          L1';
              (gamma*X1)       (-gamma*X2)  (zeros(n,n))  (zeros(n,m));
                   X1          (zeros(n,n))  -(inv(Qc))    (zeros(n,m));
                   L1          (zeros(m,n)) (zeros(m,n))    -(inv(Rc))];

LMIs=[LMIs, S_01<=0];

LMIs=[LMIs, [W1 G1'; G1 X1]>=0];
               
%% Subsistema 02

Euclideana_2=(Ah{2}*X2+Bh{2}*L2)'+(Ah{2}*X2+Bh{2}*L2);

S_02=[ (Euclideana_2-gamma*X2) (gamma*X2)'       X2'          L2';
              (gamma*X2)       (-gamma*X1)  (zeros(n,n))  (zeros(n,m));
                   X2          (zeros(n,n))  -(inv(Qc))    (zeros(n,m));
                   L2          (zeros(m,n)) (zeros(m,n))    -(inv(Rc))];  
               
LMIs=[LMIs, S_02<=0];

LMIs=[LMIs, [W2 G2'; G2 X2]>=0];
%%
options=sdpsettings;
% options.solver='sedumi';
% options.solver='sdp3';
options.solver='lmilab';
options.tol=1e-2;
options.verbose=0;
obj=max(trace(W1), trace(W2));
solvesdp(LMIs, obj, options)

X1=value(X1); X2=value(X2);
L1=value(L1); L2=value(L2);
P1=inv(X1); P2=inv(X2);

K1=L1*inv(X1);
K2=L2*inv(X2);

K1_K2=[K1', K2']
%% Fim do controlador
npts=300; ts=0.01; %npts=75;
t=linspace(0,20,2*npts);

ref=[1*ones(1,npts) 1*ones(1,npts);
     0*ones(1,npts) 2*ones(1,npts)];
 
Bz=[zeros(2,2); zeros(2,2); zeros(2,2); eye(2)];

% Entrada externa impulsiva
W=[1 zeros(1,npts-1) 1 zeros(1,npts-1);
   1 zeros(1,npts-1) 1 zeros(1,npts-1)];

%% condição inicial
x=[0 0 0 0 0 0 0 0]';

for i=1:length(t)
    
    Inc(1,i)=2*delta*rand+(1-delta);
    DELTA=eye(2)*Inc(1,i);
    
    % Regra de comutação
    sig1=x(:,i)'*P1*x(:,i);
    sig2=x(:,i)'*P2*x(:,i);
    sig=min(sig1, sig2);
    
    SIG_1(:,i)=sig1;
    SIG_2(:,i)=sig2;
    SIG(:,i)=sig;

   if sig==sig1
       A=Ah{1}; B=Bhn*DELTA; C=Ch{1}; D=Dh{1}; K=K1; G=G1; %G=rand(6,2);
   else
       A=Ah{2}; B=Bhn*DELTA; C=Ch{2}; D=Dh{2}; K=K2; G=G2; %G=rand(6,2);
   end

   u(:,i)=K*x(:,i);     
   
   % metodo de RK4
   k1=A*x(:,i) +B*u(:,i)+G*W(:,i)+Bz*ref(:,i);
   k2=A*(x(:,i) +0.5*ts*k1) +B*u(:,i)+G*W(:,i)+Bz*ref(:,i);
   k3=A*(x(:,i) +0.5*ts*k2) +B*u(:,i)+G*W(:,i)+Bz*ref(:,i);
   k4=A*(x(:,i) +ts*k3) +B*u(:,i)+G*W(:,i)+Bz*ref(:,i);
   x(:,i+1)=x(:,i)+(ts/6)*(k1+2*(k2+k3)+k4);  
   y(:,i)=C*x(:,i);
end

ISE=trace((ref-y)'*(ref-y));

fprintf('ISE = %f\n', ISE)

%% plot's
figure(3)
plot(t, y,'r','linewidth', 2); hold on
plot(t, ref,'k--','linewidth', 2); grid on
title('Resposta do sistema'); axis([0 t(end) -0.5 2.5]);
%%
figure(4);
subplot(211);
plot(t, Inc,'b','linewidth', 1); 
xlabel('time'); ylabel('\mu \in [(1-\delta) (1+\delta)]'); grid on
title('Incertezas');

subplot(212);
plot(t, SIG,'k','linewidth', 2); 
xlabel('time'); ylabel('min(\sigma_{1}, \sigma_{2})'); grid on
title('Regra de chaveamento');