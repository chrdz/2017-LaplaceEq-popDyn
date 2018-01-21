clear;

%=================== Principal eigenvalue ====================%

% Parameters :
mu=1; % choice of the function \mu

T=10^1; % for sinus function

plus=1; % choose higher value if \mu is very negetive
%(--> convergence of fixed point)

diff=0.1; % diffusion coefficient
% Be carefull not to choose the diffusion coefficient too small,
% if so, the methods solves "(\mu(x)-u)u=0", and lambda_1 seems 
% to be "-max(\mu)"

% Mesh :
N=10^3; % number of interior points --> N+2 points 
L_1=pi/T*4; % length of interval
h=1/(N+1); % space step

X=linspace(0,L_1,N+2); % space mesh

% Various functions \mu(x) (the derivative of f(x,.)) :

% piecewise constant functions :
p=10; 
% "trapz" is an Octave built-in function
% the various disp are meant to test if the functions \mu 
% have same integral

% constant :
r1=p*ones(N+2,1); 
disp(['1 : ',num2str(trapz(X,r1')),' vs ',num2str(p*L_1), ' p is ',num2str(p)])

% 2 hills :
x_0=L_1/8; x_1=L_1/8*3; x_2=L_1/8*5; x_3=L_1/8*7;
q=p*L_1/(x_1-x_0+x_3-x_2);
r2=q*(X>=x_0)-q*(X>=x_1)+q*(X>=x_2)-q*(X>=x_3);
disp(['2 : ',num2str(trapz(X,r2)),' vs ',num2str(q*4*L_1/8), ' p is ',num2str(p)])

% 5 hills :
x_0=L_1/10; x_1=L_1/10*2; x_2=L_1/10*3; x_3=L_1/10*4; 
x_4=L_1/10*5; x_5=L_1/10*6; x_6=L_1/10*7; x_7=L_1/10*8;
x_8=L_1/10*9; x_9=L_1/10*10;
q=L_1*p/(5*(x_1-x_0));
r3= q*(X>=x_0)-q*(X>=x_1)+q*(X>=x_2)-q*(X>=x_3)+q*(X>=x_4)-q*(X>=x_5)+q*(X>=x_6)-q*(X>=x_7)+q*(X>=x_8)-q*(X>=x_9);
disp(['3 : ',num2str(trapz(X,r3)),' vs ',num2str(q*5*L_1/10), ' p is ',num2str(p)])

% 2 hills + negative parts :
j=5;
x_0=L_1/8; x_1=L_1/8*3; x_2=L_1/8*5; x_3=L_1/8*7;
q=2*p+j;
r4=-j*ones(1,N+2); 
r4(X>=x_0 & X<=x_1)=q;
r4(X>=x_2 & X<=x_3)=q;
%r4=q*(X>=x_0)-q*(X>=x_1)+q*(X>=x_2)-q*(X>=x_3);
disp(['4 : ',num2str(trapz(X,r4)),' vs ',num2str(q*4*L_1/8-j*4*L_1/8), ' p is ',num2str(p)])

% 1 hill :
x_0=L_1/8*3; x_1=L_1/8*5; q=L_1*p/(x_1-x_0);
r5=q*(X>=x_0)-q*(X>=x_1);
disp(['5 : ',num2str(trapz(X,r5)),' vs ',num2str(q*2*L_1/8), ' p is ',num2str(p)])

% unfavourable zone in the middle:
y=20;k=L_1/y;
x_0=0*k; x_1=k*9; x_2=k*11; x_3=k*20;
q=p*L_1/(x_1-x_0+x_3-x_2);
r8=q*(X>=x_0)-q*(X>=x_1)+q*(X>=x_2)-q*(X>=x_3);
disp(['8 : ',num2str(trapz(X,r8)),' vs ',num2str(q*8*L_1/10), ' p is ',num2str(p)])

% Compare sinus and constant function :
r6=p*sin(T*X-L_1)+p/2; % 2 hills
fun1 = @(x) p*sin(T*x-L_1)+p/2;
int_r6 = quad(fun1,0,L_1); % computes integral
% "quad" is an Octave built-in function
p=int_r6/L_1; r7 = p*ones(N+2,1); % constant function
% tests if constant function has same integral as sine :
fun2= @(x) p; quad(fun2,0,L_1); 
disp(['6 : (sine) ',num2str(int_r6),' vs ',num2str(quad(fun2,0,L_1)), ' p is ',num2str(p)])
disp(['7 : (trapz - sine) ',num2str(trapz(X,r6)),' vs ',num2str(trapz(X,r7')), ' p is ',num2str(p)])

% Choice of the function \mu(x) :
r=eval(['r',num2str(mu)]);
c=max(r)+2;

% Discrete Laplacian :
A=sparse(1:N+2,1:N+2,2)+sparse(2:N+2,1:N+1,-1,N+2,N+2)+sparse(1:N+1,2:N+2,-1,N+2,N+2);
A(1,N+2)=-1; A(N+2,1)=-1;

% Matrix diff*(1/h^2)*A + c*I - r*I :
B=diff*1/h^2*A + sparse(1:N+2,1:N+2,c) - sparse(1:N+2,1:N+2,r);

% Computation of principal eigenvalue of the inverse of B :
lambda_1=0;
mu_1=power_iteration(B);
lambda_1=1/mu_1 - c

%=================== Stationary solution =====================%

u_k=100*rand(N+2,1); u_k1=u_k; err=4; % initialization

% I + diff*(1/h^2)*A :
C=plus*sparse(1:N+2,1:N+2,1) + diff*(1/h^2)*A; 

% fixed point method:
while (err>10^(-9))
    u_k=u_k1; % save old value
 
    D=C+sparse(1:N+2,1:N+2,u_k); % matrix we have to invert
    
    % solves D*u_k1 = I*r*u_k + u_k :
    u_k1=D\(sparse(1:N+2,1:N+2,r)*u_k+plus*u_k); 
    
    err=norm(u_k-u_k1); % computes error
end

%=================== ==== Display ============================%

% Plot of \mu :
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 (2*scrsz(3)/2)/1.2 (2*scrsz(4)/2)/2]) ; 

subplot(1,2,1); plot(X, r, 'LineWidth',2); 
set(gca,'FontSize',17);
xlim([0 X(N+2)]); ylim([-c c]);
title(['\lambda_1 = ',num2str(lambda_1)]);
xlabel('x'); ylabel('\mu(x)');

% Plot of u :
subplot(1,2,2); plot(X, u_k1,'r', 'LineWidth',2); 
set(gca,'FontSize',17); 
xlim([0 X(N+2)]); xlabel('x'); ylabel('u(x)');

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,[num2str(mu),'.png']); % saves the current figure