clc;clear;close all;
A = [ -6 -1 1; 3 -5 3;  4 1 -2];  %%aa
B = [ -2 -1 -5; 4 2 7; -2 -3 3];  %aa
C = eye(2,3);
D = zeros(2,3);
poles = [.66 .52 0.8];

% A = [ -2];  %%aa
% B = [ 3];  %aa
% C = .1*eye(1);
% D = zeros(1);
% poles = [.66];

% C = eye(3,3);
% D = zeros(3,3);
sysCont = ss(A,B,C,D);
%Discrete model
Ts = .1;
sysDiscr = c2d(sysCont,Ts,'zoh');
sys.a = sysDiscr.a;sys.b = sysDiscr.b;sys.c = sysDiscr.c;sys.d = sysDiscr.d;sys.Ts = Ts;



L_gain = place(sys.a', sys.c', poles).';

[sys.n_x,sys.n_u] = size(sys.b);
[sys.n_y,~] = size(sys.c);
sys.ap = [sys.a,sys.b;zeros(sys.n_u,sys.n_x),eye(sys.n_u)];
sys.bp = [sys.b;eye(sys.n_u)];
sys.cp = [sys.c,zeros(sys.n_y,sys.n_u)];

sys.h = 10;
sys.Q_x = 0.0*eye(sys.n_x);sys.R_u = 0.0*eye(sys.n_u);sys.R_du = 0.0*eye(sys.n_u);sys.Q_y = .4*eye(sys.n_y);
% sys.Q_x = 0.1*eye(sys.n_x);sys.R_u = 0.2*eye(sys.n_u);sys.R_du = 0.3*eye(sys.n_u);sys.Q_y = .4*eye(sys.n_y);
sys.Q_xp = blkdiag(sys.Q_x,sys.R_u);

sys.lb_u = [-100; -110; -120];
sys.ub_u = [150; 140; 130];
sys.lb_x = -Inf(sys.n_x,1);
sys.ub_x = +Inf(sys.n_x,1);
sys.lb_y = -Inf(sys.n_y,1);
sys.ub_y = +Inf(sys.n_y,1);
x_0 = zeros(sys.n_x,1);
u__1 = zeros(sys.n_u,1);
xsp = 5*ones(sys.n_x,1);
ysp = sys.c*xsp
y_sp = [];x_sp = [];
for i=1:sys.h
   y_sp = [y_sp,ysp];
   x_sp = [x_sp,xsp];
end
lb = [];ub = [];
for i=1:sys.h
    lb = [lb;sys.lb_u-sys.ub_u;sys.lb_x-x_sp(:,i);sys.lb_u];
    ub = [ub;sys.ub_u-sys.lb_u;sys.ub_x-x_sp(:,i);sys.ub_u];
end
for i=1:sys.h
    lb = [lb;sys.lb_y - y_sp(:,i)];
    ub = [ub;sys.ub_y - y_sp(:,i)];
end
x_spp = [x_sp;zeros(sys.n_u,sys.h)];
[Aeq,H,f] = constant_term(sys);
A = [];
b = [];

u__1 = zeros(sys.n_u,1);
x1_updated(:,1) = zeros(sys.n_x,1);
x_plant(:,1) = zeros(sys.n_x,1);
y1(:,1) = sys.c*x1_updated(:,1);
y_plant(:,1) = sys.c*x_plant(:,1);
NTs = 150;
umpc = zeros(size(H,1));
for i=1:NTs
    beq = b_eq_calc(sys,x_spp,x1_updated(:,i),u__1,y_sp);
    umpc = quadprog(H,f,A,b,Aeq,beq,lb,ub,umpc);                                                                %Run QP for optimal input
    u(:,i) = umpc(sys.n_u+sys.n_x+1:sys.n_u+sys.n_x+sys.n_u);                                                   %Set the first input as applied input
    u__1 = u(:,i);
    [y_plant(:,i+1),x_plant(:,i+1)] = lin_sys(x_plant(:,i),sys,u(:,i));                                         %Apply the input to the plant
    [y1(:,i+1),x1_updated(:,i+1)] = lin_sys_Luenberger_observer(x1_updated(:,i),sys,u(:,i),y_plant(:,i),L_gain);%Use plant output for state estimation    
end
T = sys.Ts*(1:NTs+1);
figure
for i=1:sys.n_y
    subplot(sys.n_y,1,i)
    %    plot(T,y_plant(i,:),[T(1),T(end)],[ysp(i),ysp(i)],'--r',T,y1(i,:))
    plot(T,y_plant(i,:),[T(1),T(end)],[ysp(i),ysp(i)],'--r')    
    legend(['y_',num2str(i)],['y_s_p_,_',num2str(i)],['y_o_b_s_v_',num2str(i)])
    xlabel('Time');ylabel(['y_',num2str(i)])
    grid on
end
figure
for i=1:sys.n_u
   subplot(sys.n_u,1,i)
   plot(T(1:end-1),u(i,:),[T(1),T(end)],[sys.lb_u(i),sys.lb_u(i)],'--r',[T(1),T(end)],[sys.ub_u(i),sys.ub_u(i)],'--r')
   legend(['u_',num2str(i)],['lb_u_',num2str(i)],['ub_u_',num2str(i)])
   xlabel('Time');ylabel(['u_',num2str(i)])
   grid on
end