close all, clear all, clc

% Article - Pradhan 2015

% Example 2
A=[ -0.764  0.387 -12.9 0 0.952  6.05;
     0.024 -0.174  4.31 0 -1.76 -0.416;
     0.006 -0.999 -0.0578 0.0369 0.0092 -0.0012;
     1 0 0 0 0 0;
     0 0 0 0 -10 0;
     0 0 0 0 0 -5];
B=[0 0 0 0 20 0;0 0 0 0 0 10]'; 
C=[0 0 1 0 0 0;0 0 0 1 0 0]; D=zeros(2);
[n,p]=size(B);
S=[C zeros(p);C*A zeros(p);zeros(p,n) eye(p)];
Qc=diag([1 1 1e5 1e5 1 1 1e6 1e6]); Rc=eye(2);
% Qc=diag([1 1 1 1 1 1 1e6 1e6]); Rc=eye(2);
z0=[0 0 0 0 0 0 1 1]';

Ah=[A zeros(size(B));-C zeros(size(D))];
Bh=[B; zeros(size(D))];
delta=0.6;
% delta = 0.9;
Delta1=diag([1+delta 1+delta]);
Delta2=diag([1-delta 1-delta]);


[Klqr,Plqr]=lqr(Ah,Bh,Qc,Rc); J=z0'*Plqr*z0,

gamma=16e5, %example 2
% gamma=14e5, %example 2

% % LMI Optimisation
opts=sdpsettings; 
opts.solver='lmilab'; opts.tol=1e-2; % tolerance 
% opts.solver='sdpt3'; opts.tol=1e-5; % tolerance
% opts.solver='sedumi'; opts.tol=1e-12; % tolerance

opts.lmilab.maxiter = 100;

opts.verbose=1; % 0 - iteraction set is omitted
                % 1 - iteraction set is apparent
               

[m,n]=size(Bh);
P=sdpvar(m,m);
Y=sdpvar(n,m);
X=[S*P; Y];
[mx,nx]=size(X);

W=sdpvar(mx,mx);
Z=sdpvar(nx,nx);

% Processo de restrição
LMIs=[P>=0];
FI=+(Ah*P+Bh*Y)+(Ah*P+Bh*Y)';
Riccati=[[FI     P          Y'; 
          P  -inv(Qc)  zeros(m,n);
          Y zeros(n,m)  -inv(Rc)  ] <= 0];
Restr5=[[gamma z0';z0 P]>=0]; 
LMIs=[LMIs,Riccati,Restr5];
obj=rank(X); % Optimal PID solution

% Suboptimal PID Solution
Restr8=[[W X;X' Z]>=0];
LMIs=[LMIs,Restr8];
obj=trace(W)+trace(Z); %Proposition 3

solvesdp(LMIs,obj,opts)
% checkset(LMIs)

P=value(P);
Y=value(Y);
K=Y*inv(P);
%%
% Kpdi=K*S'*inv(S*S');
Kpdi=K*pinv(S);

% Valores LMIs
K_p=Kpdi(:,1:2); K_d=Kpdi(:,3:4); K_i=Kpdi(:,5:6);
Kd=K_d*(eye(2)+C*B*K_d),
Kp=(eye(2)-Kd*C*B)*K_p,
Ki=(eye(2)-Kd*C*B)*K_i,

% Closed-loop time response
npts=75;
t=linspace(0,20,2*npts);
ref=[1*ones(1,npts) 1*ones(1,npts);
     0*ones(1,npts) 2*ones(1,npts)];

tau=0.01*norm(Kd)/norm(Kp);
Fs=tf(1,[tau 1],'io',30e-3);
% delta=0.6;
% Delta1=diag([1+delta 1+delta]);
% Delta2=diag([1-delta 1-delta]);

% Trace Minimized suboptimal LQR 
Bcl=[zeros(6,2); eye(2)];
Ccl=[C zeros(2)];
Dcl=zeros(2);
Tsscl=series(Fs,ss(Ah+Bh*K,Bcl,Ccl,Dcl)); [yss,t,xss]=lsim(Tsscl,ref,t);
Tsscl1=series(Fs,ss(Ah+Bh*Delta1*K,Bcl,Ccl,Dcl)); [yss1,t,xss1]=lsim(Tsscl1,ref,t);
Tsscl2=series(Fs,ss(Ah+Bh*Delta2*K,Bcl,Ccl,Dcl)); [yss2,t,xss2]=lsim(Tsscl2,ref,t);

%Indirect pseudo?inverse PID
Kpid=inv(eye(2)-Kd*C*B)*[Kp Kd Ki]*S;
Bcl=[zeros(6,2); eye(2)];
Ccl=[C zeros(2)];
Dcl=zeros(2);
Tpidcl=series(Fs,ss(Ah+Bh*Kpid,Bcl,Ccl,Dcl)); [ypid,t,xpid]=lsim(Tpidcl,ref,t);
Tpidcl1=series(Fs,ss(Ah+Bh*Delta1*Kpid,Bcl,Ccl,Dcl)); [ypid1,t,xpid1]=lsim(Tpidcl1,ref,t);
Tpidcl2=series(Fs,ss(Ah+Bh*Delta2*Kpid,Bcl,Ccl,Dcl)); [ypid2,t,xpid2]=lsim(Tpidcl2,ref,t);



figure(1)
plot(t,yss,'r',t,ypid,'b--','linewidth',3)
% legend('Trace minimized suboptimal LQR','Indirect pseudo - inverse PID')
set(gca,'fontname','Times New Roman','fontsize',12)
xlabel('Time (sec)'), ylabel('Output (volt)')
grid

figure(2)
plot(t,yss1,'r',t,ypid1,'b-.','linewidth',3)
% legend('Trace minimized suboptimal LQR','Indirect pseudo - inverse PID')
hold on
plot(t,yss2,'r',t,ypid2,'b-.','linewidth',3)

set(gca,'fontname','Times New Roman','fontsize',12)
xlabel('Time (sec)'), ylabel('Output (volt)')
grid

figure(3)
plot(t,yss+0.01*randn(size(yss)),'r',t,ypid+0.01*randn(size(yss)),'b--','linewidth',1)
% legend('Trace minimized suboptimal LQR','Indirect pseudo - inverse PID')
set(gca,'fontname','Times New Roman','fontsize',12)
xlabel('Time (sec)'), ylabel('Output (volt)')
grid

ISEpid=sum(sum((ref'-ypid)'*(ref'-ypid)));
ISEpid1=sum(sum((ref'-ypid1)'*(ref'-ypid1)));
ISEpid2=sum(sum((ref'-ypid2)'*(ref'-ypid2)));

disp('[ISE ISE1 ISE2]')
ISE=[ISEpid ISEpid1 ISEpid2]