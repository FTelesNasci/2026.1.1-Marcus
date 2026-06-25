close all
clear all
clc

%% Article - johansson 2000 => The Quadruple Tank Process
%% Example 2

A=[ -0.764  0.387 -12.9 0 0.952  6.05;
     0.024 -0.174  4.31 0 -1.76 -0.416;
     0.006 -0.999 -0.0578 0.0369 0.0092 -0.0012;
     1 0 0 0 0 0; 0 0 0 0 -10 0; 0 0 0 0 0 -5];
B=[0 0 0 0 20 0;0 0 0 0 0 10]'; 
C=[0 0 1 0 0 0;0 0 0 1 0 0]; D=zeros(2);

delta=0.60; Qc=diag([1 1 1 1 1 1 1e4 1e4]); Rc=eye(2);  gamma=1; mu=16e6;

% delta=0.85; Qc=diag([1 1 1 1 1 1 1e3 1e3]); Rc=eye(2);  gamma=1; mu=16e10;  % limite 


% delta=0.85; Qc=diag([1 1 1 1 1 1 1e3 1e3]); Rc=eye(2);  gamma=1; mu=16e6;  % limite 

Delta1=diag([1+delta 1+delta]); 
Delta2=diag([1-delta 1-delta]); 
%% Modelo aumentado
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

[n, p]=size(B); % eq.:(15)
S=[C zeros(p);C*A zeros(p);zeros(p,n) eye(p)];

%% Condições de implementação
npts=150; % npts=75;
t=linspace(0,20,2*npts);
ref=[1*ones(1,npts) 1*ones(1,npts);
     0*ones(1,npts) 2*ones(1,npts)];
Bz=[zeros(2,2); zeros(2,2); zeros(2,2); eye(2)];

[n m]=size(Bh{1});  % Tamanho da planta 

%% Parâmetros do controlador
z0=[0 0 0 0 0 0 1 1]'; %gamma=1; 

%  G1=[Delta1; Delta1; Delta1; Delta1];
%  G2=[Delta2; Delta2; Delta2; Delta2];
G1=[eye(2); eye(2); eye(2); eye(2)];
G2=[eye(2); eye(2); eye(2); eye(2)];

X1=sdpvar(n, n, 'symmetric');
X2=sdpvar(n, n, 'symmetric');
W1=sdpvar(m, m, 'symmetric');
W2=sdpvar(m, m, 'symmetric');
L1=sdpvar(m, n, 'full');
L2=sdpvar(m, n, 'full');

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

%% Modelo do Pradhan
X_1=[S*X1; L1]; X_2=[S*X2; L2];

[mx, nx]=size(X_1);

W_1=sdpvar(mx,mx); W_2=sdpvar(mx,mx);
Z_1=sdpvar(nx,nx); Z_2=sdpvar(nx,nx);

% Suboptimal PID Solution
LMIs=[LMIs, [W_1 X_1;X_1' Z_1]>=0];
LMIs=[LMIs, [W_2 X_2;X_2' Z_2]>=0];

LMIs=[LMIs, [mu z0';z0 X1]>=0]; 
LMIs=[LMIs, [mu z0';z0 X2]>=0]; 

obj_1=trace(W_1)+trace(Z_1); %Proposition 3
obj_2=trace(W_2)+trace(Z_2); %Proposition 3


obj=max(obj_1, obj_2);
% obj=max(trace(W1), trace(W2));
%% Resolução das LMIs
options=sdpsettings;
options.solver='lmilab';
% options.solver='sedumi';
% options.solver='sdpt3';

options.tol=1e-5;
options.verbose=0;
solvesdp(LMIs, obj, options)

X1=value(X1); X2=value(X2);
L1=value(L1); L2=value(L2);
P1=inv(X1); P2=inv(X2);

K1=L1*inv(X1);
K2=L2*inv(X2);

%% fim do modelo LQR Chaveado
K1', K2',

%% Transformação para PID chaveado
% Modelo para K PID 1
Kpdi1=K1*pinv(S);
K_p1=Kpdi1(:,1:2); K_d1=Kpdi1(:,3:4); K_i1=Kpdi1(:,5:6);
Kd1=K_d1*(eye(2)+C*B*Delta1*K_d1);
Kp1=(eye(2)-Kd1*C*B*Delta1)*K_p1;
Ki1=(eye(2)-Kd1*C*B*Delta1)*K_i1;

% Modelo para K PID 2
Kpdi2=K2*pinv(S);
K_p2=Kpdi2(:,1:2); K_d2=Kpdi1(:,3:4); K_i2=Kpdi2(:,5:6);
Kd2=K_d2*(eye(2)+C*B*Delta2*K_d2);
Kp2=(eye(2)-Kd2*C*B*Delta2)*K_p2;
Ki2=(eye(2)-Kd2*C*B*Delta2)*K_i2;

%% Ganhos da transformação LQR-chaveado para PID
Kpid1=inv(eye(2)-Kd1*C*B*Delta1)*[Kp1 Kd1 Ki1]*S;
Kpid2=inv(eye(2)-Kd2*C*B*Delta2)*[Kp2 Kd2 Ki2]*S;

tau=0.01*norm((Kd1+Kd2)/2)/norm((Kp1+Kp2)/2);
Fs=tf(1,[tau 1],'io',50e-3);

tau1=0.01*norm(Kd1)/norm(Kp1);
Fs1=tf(1,[tau1 1],'io',50e-3);

tau2=0.01*norm(Kd2)/norm(Kp2);
Fs2=tf(1,[tau2 1],'io',50e-3);

% Resposta no tempo do controle proposto
Tpidcl=series(Fs, ss(Ahn+Bhn*0.5*(Kpid1+Kpid2), Bz, Chn, Dhn)); [ypid,  t, xpid]=lsim(Tpidcl,ref,t);
Tpidcl1=series(Fs1,ss(Ahn+Bhn*Delta1*Kpid1, Bz, Chn, Dhn));     [ypid1, t, xpid1]=lsim(Tpidcl1,ref,t);
Tpidcl2=series(Fs2,ss(Ahn+Bhn*Delta2*Kpid2, Bz, Chn, Dhn));     [ypid2, t, xpid2]=lsim(Tpidcl2,ref,t);
%% Ganho do PID proposto por Pradhan 2015

%% Análise comparativa com o modelo do pradhan 2015
Kp_p=[-5.030 -0.253;2.830 -3.457];
Ki_p=[ 6.667 0.451;-0.995 4.098];
Kd_p=[-2.183 -0.490;0.758 -1.281];

tau=0.01*norm(Kd_p)/norm(Kp_p);
Fs=tf(1,[tau 1],'io',30e-3);

Kpidp=inv(eye(2)-Kd_p*C*B)*[Kp_p Kd_p Ki_p]*S;

% Resposta no tempo Pradhan
Tpidclp=series(Fs,ss(Ahn+Bhn*Kpidp, Bz, Chn, Dhn));         [ypidp,t,xpid]=lsim(Tpidclp,ref,t);
Tpidcl1p=series(Fs,ss(Ahn+Bhn*Delta1*Kpidp, Bz, Chn, Dhn)); [ypid1p,t,xpid1p]=lsim(Tpidcl1p,ref,t);
Tpidcl2p=series(Fs,ss(Ahn+Bhn*Delta2*Kpidp, Bz, Chn, Dhn)); [ypid2p,t,xpid2p]=lsim(Tpidcl2p,ref,t);

%% indices de desempenho
ISEpid=trace((ref'-ypid)'*(ref'-ypid));
ISEpid1=trace((ref'-ypid1)'*(ref'-ypid1));
ISEpid2=trace((ref'-ypid2)'*(ref'-ypid2));

ISEpidp=trace((ref'-ypidp)'*(ref'-ypidp));
ISEpid1p=trace((ref'-ypid1p)'*(ref'-ypid1p));
ISEpid2p=trace((ref'-ypid2p)'*(ref'-ypid2p));

disp('Proposto -> [ISE ISE1 ISE2]')
ISE=[ISEpid ISEpid1 ISEpid2]
disp('Pradhan -> [ISEp ISE1p ISE2p]')
ISEp=[ISEpidp ISEpid1p ISEpid2p]

%% Plots 01

figure(1)
subplot(2,1,1)
plot(t,ypid1(:,1),'r',t,ypid1p(:,1),'b',t,ref(1,:)','k--',...
     t,ypid1(:,2),'r',t,ypid1p(:,2),'b',t,ref(2,:)','k--','linewidth',1.5), 
axis([0 t(end) -.5 2.5]), grid; title('Ponto de operação 01 (1+\delta)')
legend('Proposed','Benchmark','Setpoint'); xlabel('Time(sec)'); ylabel('Output(rad)');

subplot(2,1,2)
plot(t,ypid2(:,1),'r',t,ypid2p(:,1),'b',t,ref(1,:)','k--',...
     t,ypid2(:,2),'r',t,ypid2p(:,2),'b',t,ref(2,:)','k--','linewidth',1.5), 
axis([0 t(end) -1 3]); grid; title('Ponto de operação 02 (1-\delta)'); 
% legend('Proposed','Benchmark','Setpoint')
 xlabel('Time(sec)'); ylabel('Output(rad)');