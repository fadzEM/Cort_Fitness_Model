
close all
%s = sym('s');
%parameters
beta = 1; % hill function parameter
K = 5; % hill function half saturation constant
T = 2;% age at which birds start breeding
b3 = 2.5; % decay rate for group 3
b1 = 3; % decay rate for group 1
b2 = 1.5; % decay rate for group 2
B = 15; % constant growth rate for bird population is alpha*B
alpha1 = .2; % alpha for group 1
alpha3 = .5; % alpha for group 3
alpha2 = .8; % alpha for group 2
c_T = .2; % Cort threshold
c_0 = .05; % Cort minimal value
m0 = .6;% birds natural mortality rate


%  cort vs environmental challenges
N = 1;
x = 0:.01:N; %Environmental Challenges

yTrue = zeros(N,length(x));
yTrue(6,:) = cort(beta, K,6,x);
yTrue(1,:) = cort(beta, K,1,x);
yTrue(10,:) = cort(beta, K,10,x);
for i = 1:11
    yTrue(i+2,:) = cort(beta, K,i+2,x);
end    

% Cort vs EC figure cortvsEC
% figure(1)
% plot(x,yTrue(1,:), 'k:', 'LineWidth', 2)
% hold on
% plot(x,yTrue(3,:), 'k--', 'LineWidth', 2)
% hold on
% plot(x,yTrue(9,:), 'k', 'LineWidth', 2)
% xlim([0, 1.1])
% ylim([0, 1.1])
% xlabel('Environmental Challenges', 'Fontsize', 20)
% ylabel('Cort Level', 'Fontsize', 20)
% h=legend( 'predictable', 'semi-predictable', 'unpredictable');%,'Fontsize', 20)%, 'n=6','n=7','n=8','n=9','n=10')
% set(h,'FontSize',20);


% reproductive success vs cort for all groups 
c = 0:0.01:1;
F1 = bx(alpha1,B,b1,c, c_T);
F2 = bx(alpha2,B,b2,c, c_T);
F3 = bx(alpha3,B,b3,c, c_T);
figure(2)
plot(c,F1/max(F2), 'k:', 'LineWidth', 2)
hold on
plot(c,F2/max(F2), 'k--', 'LineWidth', 2)
hold on
plot(c,F3/max(F2), 'k', 'LineWidth', 2)
hold on
plot(c_T,F1/max(F2), 'r', 'LineWidth', 2)
ylim([0, 1.1])
xlabel('Cort Level', 'Fontsize', 20)
ylabel('Reproductive Success (b_x/max(b_x))', 'Fontsize', 20)
h=legend( 'Group 1', 'Group 2', 'Group 3');%,'Fontsize', 20)%, 'n=6','n=7','n=8','n=9','n=10')
set(h,'FontSize',20);
plottools


% reproductive success vs EC (unpredictalbe) for all groups
E1 = bx(alpha1,B,b1,yTrue(9,:), c_T);
E2 = bx(alpha2,B,b2,yTrue(9,:), c_T);
E3 = bx(alpha3,B,b3,yTrue(9,:), c_T);
% figure(3)
% plot(x,E1/max(E2), 'k:', 'LineWidth', 1.5)
% hold on
% plot(x,E2/max(E2), 'k--', 'LineWidth', 1.5)
% hold on
% plot(x,E3/max(E2), 'k', 'LineWidth', 1.5)
% ylim([0, 1.1])
% xlabel('Environmental Challenges', 'Fontsize', 20)
% ylabel('Reproductive Success (b_x/max(b_x))', 'Fontsize', 20)
% h=legend( 'Group 1', 'Group 2', 'Group 3');%,'Fontsize', 20)%, 'n=6','n=7','n=8','n=9','n=10')
% set(h,'FontSize',20);
% plottools


% Solving for r using Newton method.
f = @(T,bxx,m0) log(abs(s+m0))-log(bxx) + (s+m0).*T;
df = @(T,m0) T + (1./abs(s+m0));


R1 =zeros(1,length(x));
R2 = zeros(1,length(x));
R3 = zeros(1,length(x));

K1 =zeros(1,length(x));
K2 = zeros(1,length(x));
K3 = zeros(1,length(x));

for j = 1:length(x)
    [R1all, exR1] = newton(f(2,F1(j),m0),df(T,m0),.1, 0.5*10^-5, 10);R1(j) = R1all(end);
    [R2all, exR2] = newton(f(2,F2(j),m0),df(T,m0),.1, 0.5*10^-5, 10);R2(j) = R2all(end);
    [R3all, exR3] = newton(f(2,F3(j),m0),df(T,m0),0.1, 0.5*10^-5, 10);R3(j) = R3all(end);
    [K1all, exK1] = newton(f(2,E1(j),m0),df(T,m0),.1, 0.5*10^-5, 10);K1(j) = K1all(end);
    [K2all, exK2] = newton(f(2,E2(j),m0),df(T,m0),.1, 0.5*10^-5, 10);K2(j) = K2all(end);
    [K3all, exK3] = newton(f(2,E3(j),m0),df(T,m0),0.1, 0.5*10^-5, 10);K3(j) = K3all(end);
end    

% r vs cort for all groups
figure(4)
plot(c,R1/max(R2), 'k:', 'LineWidth', 2)
hold on
plot(c,R2/max(R2), 'k--', 'LineWidth', 2)
hold on
plot(c,R3/max(R2), 'k', 'LineWidth', 2)
hold on
plot(c,0*R1, 'r', 'LineWidth', 2)
ylim([-0.8, 1.1])
xlabel('Baseline Cort', 'Fontsize', 20)
ylabel('Intrinsic Rate of Natural Increase (r/max(r))', 'Fontsize', 20)
h=legend( 'Group 1', 'Group 2', 'Group 3');%,'Fontsize', 20)%, 'n=6','n=7','n=8','n=9','n=10')
set(h,'FontSize',20);
plottools

% r vs EC (unpredictable) for all groups
figure(5)
plot(x,K1/max(K2), 'k:', 'LineWidth', 2)
hold on
plot(x,K2/max(K2), 'k--', 'LineWidth', 2)
hold on
plot(x,K3/max(K2), 'k', 'LineWidth', 2)
hold on
plot(x,0*K1,'r','LineWidth', 2)
ylim([-0.8, 1.1])
xlabel('Environmental Challenges', 'Fontsize', 20)
ylabel('Intrinsic Rate of Natural Increase (r/max(r))', 'Fontsize', 20)
h=legend( 'Group 1', 'Group 2', 'Group 3');%,'Fontsize', 20)%, 'n=6','n=7','n=8','n=9','n=10')
set(h,'FontSize',20);
plottools



% r vs EC for each group. Figure randfvsEC
% for predictable EC
bx11 = bx(alpha1,B,b1,yTrue(1,:),c_T);
bx21 = bx(alpha2,B,b2,yTrue(1,:),c_T);
bx31 = bx(alpha3,B,b3,yTrue(1,:),c_T);

r11 =zeros(1,length(x));
r21 = zeros(1,length(x));
r31 = zeros(1,length(x));
for j = 1:length(x)
    [r11all, ex11] = newton(f(2,bx11(j),m0),df(T,m0),.1, 0.5*10^-5, 10);r11(j) = r11all(end);%r(2,yTrue(1,:),c_T,m0, alpha1, b1, B);
    [r21all, ex21] = newton(f(2,bx21(j),m0),df(T,m0),.1, 0.5*10^-5, 10);r21(j) = r21all(end);
    [r31all, ex31] = newton(f(2,bx31(j),m0),df(T,m0),0.1, 0.5*10^-5, 10);r31(j) = r31all(end);
end    
f11= bx(alpha1,B,b1,yTrue(1,:), c_T);
f21= bx(alpha2,B,b2,yTrue(1,:), c_T);
f31= bx(alpha3,B,b3,yTrue(1,:), c_T);

% % for semi-predictable EC
bx13 = bx(alpha1,B,b1,yTrue(3,:),c_T);
bx23 = bx(alpha2,B,b2,yTrue(3,:),c_T);
bx33 = bx(alpha3,B,b3,yTrue(3,:),c_T);

r13 =zeros(1,length(x));
r23 = zeros(1,length(x));
r33 = zeros(1,length(x));
for j = 1:length(x)
    [r13all, ex11] = newton(f(2,bx13(j),m0),df(T,m0),.1, 0.5*10^-5, 10);r13(j) = r13all(end);%r(2,yTrue(1,:),c_T,m0, alpha1, b1, B);
    [r23all, ex21] = newton(f(2,bx23(j),m0),df(T,m0),.1, 0.5*10^-5, 10);r23(j) = r23all(end);
    [r33all, ex31] = newton(f(2,bx33(j),m0),df(T,m0),0.1,0.5*10^-5, 10);r33(j) = r33all(end);
end    

f13= bx(alpha1,B,b1,yTrue(3,:), c_T);
f23= bx(alpha2,B,b2,yTrue(3,:), c_T);
f33= bx(alpha3,B,b3,yTrue(3,:), c_T);
% 
% for Unpredictable EC
bx19 = bx(alpha1,B,b1,yTrue(9,:),c_T);
bx29 = bx(alpha2,B,b2,yTrue(9,:),c_T);
bx39 = bx(alpha3,B,b3,yTrue(9,:),c_T);

r19 =zeros(1,length(x));
r29 = zeros(1,length(x));
r39 = zeros(1,length(x));
for j = 1:length(x)
    [r19all, ex19] = newton(f(2,bx19(j),m0),df(T,m0),.1, 0.5*10^-5, 10);r19(j) = r19all(end);%r(2,yTrue(1,:),c_T,m0, alpha1, b1, B);
    [r29all, ex29] = newton(f(2,bx29(j),m0),df(T,m0),.1, 0.5*10^-5, 10);r29(j) = r29all(end);
    [r39all, ex39] = newton(f(2,bx39(j),m0),df(T,m0),0.1,0.5*10^-5, 10);r39(j) = r39all(end);
end    

f19= bx(alpha1,B,b1,yTrue(9,:), c_T);
f29= bx(alpha2,B,b2,yTrue(9,:), c_T);
f39= bx(alpha3,B,b3,yTrue(9,:), c_T);


% 
figure(6)
subplot(3,2,1)
plot(x,r11/abs(max(r21)), 'k:', 'LineWidth', 2)
hold on
plot(x,r21/abs(max(r21)), 'k--', 'LineWidth', 2)
hold on
plot(x,r31/abs(max(r21)), 'k', 'LineWidth', 2)
hold on
plot(x,0*r21,'r', 'LineWidth', 2);
ylim([-0.6, 1.1])
ylabel('r/max(r) for Predictable', 'Fontsize', 20)

subplot(3,2,3)
plot(x,r13/abs(max(r23)), 'k:', 'LineWidth', 2)
hold on
plot(x,r23/abs(max(r23)), 'k--', 'LineWidth', 2)
hold on
plot(x,r33/abs(max(r23)), 'k', 'LineWidth', 2)
hold on
plot(x,0*r13,'r', 'LineWidth', 2);
ylim([-0.8, 1.1])
ylabel('r for Semi-Predictable', 'Fontsize', 20)
subplot(3,2,5)
plot(x,r19/abs(max(r29)), 'k:', 'LineWidth', 2)
hold on
plot(x,r29/abs(max(r29)), 'k--', 'LineWidth', 2)
hold on
plot(x,r39/abs(max(r29)), 'k', 'LineWidth', 2)
hold on
plot(x,0*r19,'r', 'LineWidth', 2);
ylim([-1.1, 1.1])
xlabel('Environmental Challenges', 'Fontsize', 20)
ylabel('r for Unpredictable', 'Fontsize', 20)
h2=legend( 'Group 1','Group 2','Group 3');
set(h2,'FontSize',20);
subplot(3,2,2)
plot(x,f19/abs(max(f29)), 'k:', 'LineWidth', 2)
hold on
plot(x,f29/abs(max(f29)), 'k--', 'LineWidth', 2)
hold on
plot(x,f39/abs(max(f29)), 'k', 'LineWidth', 2)
hold on
plot(x,0*f19,'r', 'LineWidth', 2);
ylim([-0.1, 1.1])
ylabel('b_x for Predictable', 'Fontsize', 20)

subplot(3,2,4)
plot(x,f13/abs(max(f23)), 'k:', 'LineWidth', 2)
hold on
plot(x,f23/abs(max(f23)), 'k--', 'LineWidth', 2)
hold on
plot(x,f33/abs(max(f23)), 'k', 'LineWidth', 2)
hold on
plot(x,0*f13,'r', 'LineWidth', 2);
ylim([-0.1, 1.1])
ylabel('b_x for Semi-Predictable', 'Fontsize', 20)
subplot(3,2,6)
plot(x,f19/abs(max(f29)), 'k:', 'LineWidth', 2)
hold on
plot(x,f29/abs(max(f29)), 'k--', 'LineWidth', 2)
hold on
plot(x,f39/abs(max(f29)), 'k', 'LineWidth', 2)
hold on
plot(x,0*f19,'r', 'LineWidth', 2);
ylim([-0.1, 1.1])
xlabel('Environmental Challenges', 'Fontsize', 20)
ylabel('b_x for Unpredictable', 'Fontsize', 20)
plottools



% r and bx vs cort for all three groups
figure(7)
subplot(3,2,1)
plot(yTrue(1,:),r11/abs(max(r21)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(1,:),r21/abs(max(r21)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(1,:),r31/abs(max(r21)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(1,:),0*r11,'r', 'LineWidth', 2);
ylim([-0.6, 1.1]);
xlim([0.06, 0.72]);
ylabel('r for Predictable', 'Fontsize', 20)

subplot(3,2,3)
plot(yTrue(3,:),r13/abs(max(r23)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(3,:),r23/abs(max(r23)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(3,:),r33/abs(max(r23)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(3,:),0*r13,'r', 'LineWidth', 2);
ylim([-0.8, 1.1]);
xlim([0.06, 0.92]);
ylabel('r for Semi-Predictable', 'Fontsize', 20)
subplot(3,2,5)
plot(yTrue(9,:),r19/abs(max(r29)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(9,:),r29/abs(max(r29)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(9,:),r39/abs(max(r29)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(9,:),0*r19,'r', 'LineWidth', 2);
ylim([-1.1, 1.1]);
xlim([0.06, 1]);
xlabel('Cort', 'Fontsize', 20)
ylabel('r for Unpredictable', 'Fontsize', 20)
h2=legend( 'Group 1','Group 2','Group 3');
set(h2,'FontSize',20);
subplot(3,2,2)
plot(yTrue(1,:),f11/abs(max(f21)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(1,:),f21/abs(max(f21)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(1,:),f31/abs(max(f21)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(1,:),0*f19,'r', 'LineWidth', 2);
ylim([-0.1, 1.1])
xlim([0.06, 0.72]);
ylabel('b_x for Predictable', 'Fontsize', 20)

subplot(3,2,4)
plot(yTrue(3,:),f13/abs(max(f23)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(3,:),f23/abs(max(f23)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(3,:),f33/abs(max(f23)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(3,:),0*f13,'r', 'LineWidth', 2);
ylim([-0.1, 1.1])
xlim([0.06, 0.92]);
ylabel('b_x for Semi-Predictable', 'Fontsize', 20)
subplot(3,2,6)
plot(yTrue(9,:),f19/abs(max(f29)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(9,:),f29/abs(max(f29)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(9,:),f39/abs(max(f29)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(9,:),0*f19,'r', 'LineWidth', 2);
ylim([-0.1, 1.1])
xlim([0.06, 1]);
xlabel('Cort', 'Fontsize', 20)
ylabel('b_x for Unpredictable', 'Fontsize', 20)
plottools





