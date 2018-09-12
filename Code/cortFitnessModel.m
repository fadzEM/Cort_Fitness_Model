
close all

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
m0 = .09;% birds natural mortality rate


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
figure(1)
plot(x,yTrue(1,:), 'k:', 'LineWidth', 2)
hold on
plot(x,yTrue(3,:), 'k--', 'LineWidth', 2)
hold on
plot(x,yTrue(9,:), 'k', 'LineWidth', 2)
xlim([0, 1.1])
ylim([0, 1.1])
xlabel('Environmental Challenges', 'Fontsize', 20)
ylabel('Cort Level', 'Fontsize', 20)
h=legend( 'predictable', 'semi-predictable', 'unpredictable');%,'Fontsize', 20)%, 'n=6','n=7','n=8','n=9','n=10')
set(h,'FontSize',20);


% %cort effect on reproduction figure RepRate
c = 0:0.01:1;
F = f(alpha1,B,b1,c, c_T);
figure(2)
plot(c,F/max(F), 'k', 'LineWidth', 1.5)
ylim([0, 1.1])
xlabel('Baseline Cort', 'Fontsize', 20)
ylabel('Number of Offspring ', 'Fontsize', 20)

% r vs EC for B1, B2, B3 at n1, n3, n9 figure randfvsEC
% for group 1
r11 = fitness(2,yTrue(1,:),c_T,m0, alpha1, b1, B);
r21 = fitness(2,yTrue(1,:),c_T,m0, alpha2, b2, B);
r31 = fitness(2,yTrue(1,:),c_T,m0, alpha3, b3, B);
f11= f(alpha1,B,b1,yTrue(1,:), c_T);
f21= f(alpha2,B,b2,yTrue(1,:), c_T);
f31= f(alpha3,B,b3,yTrue(1,:), c_T);
% 
% % for group 3
r13 = fitness(2,yTrue(3,:),c_T,m0, alpha1, b1, B);
r23 = fitness(2,yTrue(3,:),c_T,m0, alpha2, b2, B);
r33 = fitness(2,yTrue(3,:),c_T,m0, alpha3, b3, B);
f13= f(alpha1,B,b1,yTrue(3,:), c_T);
f23= f(alpha2,B,b2,yTrue(3,:), c_T);
f33= f(alpha3,B,b3,yTrue(3,:), c_T);
% 
% % for group 2
r19 = fitness(2,yTrue(9,:),c_T,m0, alpha1, b1, B);
r29 = fitness(2,yTrue(9,:),c_T,m0, alpha2, b2, B);
r39 = fitness(2,yTrue(9,:),c_T,m0, alpha3, b3, B);
f19= f(alpha1,B,b1,yTrue(9,:), c_T);
f29= f(alpha2,B,b2,yTrue(9,:), c_T);
f39= f(alpha3,B,b3,yTrue(9,:), c_T);

% 
figure(3)
subplot(3,2,1)
plot(x,r11/abs(max(r21)), 'k:', 'LineWidth', 2)
hold on
plot(x,r21/abs(max(r21)), 'k--', 'LineWidth', 2)
hold on
plot(x,r31/abs(max(r21)), 'k', 'LineWidth', 2)
hold on
plot(x,0*r11,'r', 'LineWidth', 2);
ylim([-0.5, 1.1])
ylabel('r for Predictable', 'Fontsize', 20)

subplot(3,2,3)
plot(x,r13/abs(max(r23)), 'k:', 'LineWidth', 2)
hold on
plot(x,r23/abs(max(r23)), 'k--', 'LineWidth', 2)
hold on
plot(x,r33/abs(max(r23)), 'k', 'LineWidth', 2)
hold on
plot(x,0*r13,'r', 'LineWidth', 2);
ylim([-0.5, 1.1])
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
plot(x,0*f11,'r', 'LineWidth', 2);
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




% r and bx vs cort for all three groups
figure(4)
subplot(3,2,1)
plot(yTrue(1,:),r11/abs(max(r21)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(1,:),r21/abs(max(r21)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(1,:),r31/abs(max(r21)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(1,:),0*r11,'r', 'LineWidth', 2);
ylim([-0.5, 1.1])
ylabel('r for Predictable', 'Fontsize', 20)

subplot(3,2,3)
plot(yTrue(3,:),r13/abs(max(r23)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(3,:),r23/abs(max(r23)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(3,:),r33/abs(max(r23)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(3,:),0*r13,'r', 'LineWidth', 2);
ylim([-0.5, 1.1])
ylabel('r for Semi-Predictable', 'Fontsize', 20)
subplot(3,2,5)
plot(yTrue(9,:),r19/abs(max(r29)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(9,:),r29/abs(max(r29)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(9,:),r39/abs(max(r29)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(9,:),0*r19,'r', 'LineWidth', 2);
ylim([-1.1, 1.1])
xlabel('Cort', 'Fontsize', 20)
ylabel('r for Unpredictable', 'Fontsize', 20)
h2=legend( 'Group 1','Group 2','Group 3');
set(h2,'FontSize',20);
subplot(3,2,2)
plot(yTrue(9,:),f19/abs(max(f29)), 'k:', 'LineWidth', 2)
hold on
plot(yTrue(9,:),f29/abs(max(f29)), 'k--', 'LineWidth', 2)
hold on
plot(yTrue(9,:),f39/abs(max(f29)), 'k', 'LineWidth', 2)
hold on
plot(yTrue(9,:),0*f19,'r', 'LineWidth', 2);
ylim([-0.1, 1.1])
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
xlabel('Cort', 'Fontsize', 20)
ylabel('b_x for Unpredictable', 'Fontsize', 20)






% % r vs age with fix EC
Yn1 = zeros(3,3);
Yn3 = zeros(3,3);
Yn9 = zeros(3,3);
Fn1 = zeros(3,3);
Fn3 = zeros(3,3);
Fn9 = zeros(3,3);
EC_new = [25, 40, 50];
alpha = [alpha1, alpha2, alpha3];
b = [b1, b2, b3];
for j = 1:3
    for i = 1:3
        Yn1(j,i) = fitness(2,yTrue(1,EC_new(j)),c_T,m0, alpha(i), b(i), B); 
        Yn3(j,i) = fitness(2,yTrue(4,EC_new(j)),c_T,m0, alpha(i), b(i), B); 
        Yn9(j,i) = fitness(2,yTrue(13,EC_new(j)),c_T,m0, alpha(i), b(i), B);
        Fn1(j,i) = f(alpha(i),B,b(i),yTrue(1,EC_new(j)), c_T); 
        Fn3(j,i) = f(alpha(i),B,b(i),yTrue(2,EC_new(j)), c_T); 
        Fn9(j,i) = f(alpha(i),B,b(i),yTrue(13,EC_new(j)), c_T);
    end
end    


% r vs age
figure(5)
plot(Yn9(3,:), 'k:', 'LineWidth', 2)
hold on
plot(Yn9(3,2), 'k--', 'LineWidth', 2)
hold on
plot(Yn9(3,3), 'k', 'LineWidth', 2)
hold on
plot(0*Yn9(3,:), 'r', 'LineWidth', 2)
hold on
scatter(1,Yn9(3,1), 'k*', 'LineWidth', 2);
hold on
scatter(2,Yn9(3,2), 'k*', 'LineWidth', 2);
hold on
scatter(3,Yn9(3,3), 'k*', 'LineWidth', 2);
ylabel('r', 'Fontsize', 20)
xlim([0.9, 3.2])
ylim([-.4, 1])
h1=legend( 'Unpredicted Event');
set(h1,'FontSize',20);

