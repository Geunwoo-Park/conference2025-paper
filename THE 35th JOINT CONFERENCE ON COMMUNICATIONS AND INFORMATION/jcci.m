clc; clear; close all;
count_mat=zeros(1,11);
cost_mat=zeros(1,11);
%for i=0:10

%% ETC LMI

A = [1 0.1; 0 2]; 

k=0.787;
B = [0; 0.1*k];  
C = [1, 0];

Dc=0;

T_s = 0.1; % 샘플링 시간

%[A,B,C,D]=c2dm(Ac,Bc,Cc,Dc,T_s);

A_bar=[A B;0 0 1];
C_bar=[C 0];
x0=[0;0];

desired_poles_Kc = [0.3, 0.4];  % 제어기 극점
desired_poles_Kobar = [0.6, 0.7,0.8]; 

Kc = -place(A, B, desired_poles_Kc);

% 관측기 게인 L 계산 (A', C'을 사용하여 관측기 설계)
K_obar = place(A_bar', C_bar', desired_poles_Kobar)';

%가중치
Q=1;
R=0.0000002;

M=[1 0 0;0 1 0];
N_1=[0 0 1];

%% 시스템 구현

T = 10;   % 총 시뮬레이션 시간
N = T / T_s; % 반복 횟수

x = zeros(2, N+1); % 실제 상태
x_hat = zeros(2, N+1); % 추정 상태
x_til = zeros(2, N+1); % 상태 오차
d_hat = zeros(1, N+1); % 외란 추정값
d_til = zeros(1, N+1); % 외란 오차값
y = zeros(1, N+1); % 출력
xt_hat = zeros(2, N+1); % 트리거링 추정 상태
dt_hat = zeros(1, N+1); % 트리거링 외란 추정값

x_eq=zeros(2, N+1);
x_bar=zeros(2, N+1);
episilon=zeros(2, N+1);
sigma=zeros(1, N+1); 

x(:,1) = [0.1; 0]; % 초기 상태
x_hat(:,1) = [0; 0]; % 초기 추정 상태
x_til(:,1) = [0; 0]; % 초기 오차 상태
d_hat(1) = 0; % 초기 외란 추정
d_til(:,1) = 0; % 초기 오차 상태
y(1) = 0; % 초기 출력
um(1) =0; %초기 입력
xt_hat(:,1) = [0; 0]; % 초기 트리거링 추정 상태
dt_hat(1) = 0; % 초기 트리거링 외란 추정

x_eq(:,1)=[0; 0];
x_bar(:,1)=[0; 0];
episilon(:,1)=[0; 0];
sigma(1)=0;


% 지령치 및 외란 함수
r = @(t) 1*double(t >= 3);  % 지령치 r(t)
d = @(t)  1*double(t >= 1);  % 외란 d(t)

% 시뮬레이션 루프
for k = 1:N
   t = (k-1) * T_s;

 x_til(:,k)=x(:,k)-x_hat(:,k);
 d_til(:,k)=d(t)-d_hat(:,k);
    
%변수 설정
P_1=sdpvar(3);
P_2=sdpvar(2);

Y=sdpvar(1,2);
X=sdpvar(3,1);

Qc=  [15 0;0 15];
Qe=   [65.2961    1.8149;
    1.8149   31.6055];
Qd=  10;

%LMI
%1행
M11=(A+B*Kc)'*P_2*(A+B*Kc)-P_2+Kc'*R*Kc+Qc;
M12=(A+B*Kc)'*P_2*(B*Kc)+Kc'*R*Kc;
M13=(A+B*Kc)'*P_2*(-B)+Kc'*R;
M14=(A+B*Kc)'*P_2*(M*K_obar*C_bar)+Kc'*R*N_1;

%2행
M22=(B*Kc)'*P_2*(B*Kc)+Kc'*R*Kc-Qe;
M23=(B*Kc)'*P_2*(-B)+Kc'*R;
M24=(B*Kc)'*P_2*(M*K_obar*C_bar)+Kc'*R*N_1;

%3행
M33=(-B)'*P_2*(-B)+R-Qd;
M34=(-B)'*P_2*(M*K_obar*C_bar)+R*N;

%4행
M44=(A_bar-K_obar*C_bar)'*P_1*(A_bar-K_obar*C_bar)-P_1+(M*K_obar*C_bar)'*(M*K_obar*C_bar)+Q+N_1'*R*N_1;


sig=[M11 M12 M13 M14;
    M12' M22 M23 M24;
    M13' M23' M33 M34;
    M14' M24' M34' M44 ];


cond1 = sig <= -1e-6; % 제약 조건
cond3= P_1 >= 1e-6;
cond4= P_2 >= 1e-6;

total_cond = [cond1,cond3,cond4];

optimize(total_cond); 

P_1=double(P_1);
P_2=double(P_2);
sig=double(sig);


%gain 계산
Ko=K_obar(1:2);
Kd=K_obar(3);
Kr=inv(C*inv(eye(2)-(A+B*Kc))*B); %feedforwad gain

    e = C * (x(:,k) - x_hat(:,k));   % 오차 계산
    u = Kc * xt_hat(:,k) + Kr * r(t) - dt_hat(k); % 제어 입력
    um(k+1)=u-Kc*x_eq(:,k)-Kr*r(t)+d(t);

    % 이산시간 상태 업데이트
    x(:,k+1) = A * x(:,k) + B * (u + d(t));   % 실제 시스템
    x_hat(:,k+1) = A * x_hat(:,k) + B* (u + d_hat(k)) + Ko * e; % 상태 관측기
    d_hat(k+1) = d_hat(k) + Kd * e;  % 외란 추정
    y(k+1)=C*x(:,k+1);


    % 이벤트 트리거링
    x_eq(:,k+1)=inv(eye(2)-(A+B*Kc))*B*Kr*r(t);
    x_bar(:,k+1)=x_hat(:,k+1)-x_eq(:,k+1);
    episilon(:,k+1)=xt_hat(:,k)-x_hat(:,k+1);
    sigma(k+1)=dt_hat(k)-d_hat(k+1);

    % 트리거링 조건
    if episilon(:,k+1)'*Qe*episilon(:,k+1)+sigma(k+1)'*Qd*sigma(k+1)'-x_bar(:,k+1)'*Qc*x_bar(:,k+1) >=0
           xt_hat(:,k+1)=x_hat(:,k+1);  % 송신
           dt_hat(:,k+1)=d_hat(:,k+1);
        else 
            xt_hat(:,k+1) = xt_hat(:,k);  % 이전 값을 유지
            dt_hat(:,k+1) = dt_hat(:,k); 
    end
end

t = 0:T_s:T;
%%
figure(1)
plot(t,x(1,:),'Color',"#D95319",'LineWidth',1)
title('(a)','FontName', 'Times New Roman', 'FontSize', 20);
hold on
grid on
plot(t,x_hat(1,:),'Color',"#EDB120",'LineWidth',1)
plot(t,xt_hat(1,:),'Color',"#0072BD",'LineWidth',1)
xlabel('Time (s)', 'FontName', 'Times New Roman','FontSize', 20)
legend('$x(t)$','$\hat{x}_1(t)$','$\hat{x}_1(t_k)$','Interpreter','latex', 'FontName', 'Times New Roman','FontSize', 20)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 

figure(2)
plot(t,x(2,:),'Color',"#D95319",'LineWidth',1)
title('(b)','FontName', 'Times New Roman', 'FontSize', 20);
hold on
grid on
plot(t,x_hat(2,:),'Color',"#EDB120",'LineWidth',1)
plot(t,xt_hat(2,:),'Color',"#0072BD",'LineWidth',1)
xlabel('Time (s)', 'FontName', 'Times New Roman','FontSize', 20)
legend('$x_2(t)$','$\hat{x}_2(t)$','$\hat{x}_2(t_k)$','Interpreter','latex', 'FontName', 'Times New Roman','FontSize', 20)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 

figure(3)
plot(t,d(t),'Color',"#D95319",'LineWidth',1)
title('(c)','FontName', 'Times New Roman', 'FontSize', 20);
hold on
grid on
plot(t,d_hat(1,:),'Color',"#EDB120",'LineWidth',1)
plot(t,dt_hat(1,:),'Color',"#0072BD",'linewidth',1)
xlabel('Time (s)', 'FontName', 'Times New Roman','FontSize', 20)
legend('$d(t)$','$\hat{d}(t)$','$\hat{d}(t_k)$','Interpreter','latex', 'FontName', 'Times New Roman','FontSize', 20)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 

figure(4)
plot(t,r(t),'Color',"#0072BD",'LineWidth',1)
title('(d)','FontName', 'Times New Roman', 'FontSize', 20);
hold on
grid on
plot(t,y,'Color',"#D95319",'LineWidth',1)
xlabel('Time (s)', 'FontName', 'Times New Roman','FontSize', 20)
legend('$r(t)$','$y(t)$','Interpreter','latex', 'FontName', 'Times New Roman','FontSize', 20)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 

figure(5)
plot(t,um,'Color',"#0072BD",'LineWidth',1)
title('(e)','FontName', 'Times New Roman', 'FontSize', 20);
hold on
grid on
xlabel('Time (s)', 'FontName', 'Times New Roman','FontSize', 20)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 

%% simulink용
%out=sim("event_trigger_simulink.mdl");

%count_mat(i+1) = sum(out.sample.Data >= 0)
%x_count = 1:length(count_mat);
%plot(x_count,count_mat)

%t_cost=out.cost.time;
%data_cost=out.cost.Data;

%c_count = 1:length(count_mat);
%cost_mat(i+1)=trapz(t_cost, data_cost)
% end