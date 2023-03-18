clear all

% model parameters in 1/hours: 
% r is stem division rate, 
% l stands for lambda, TA division rate
% g stands for gamma, FD apoptosis rate 

l=1/14;
r=1/(24);
g=1/36;

% average number of cells per crypt
ntot=2392.10;

% for-loop for different number of stem cells
X1=zeros(1601,1);
P_TA=zeros(1601,5);
X2=zeros(5001,1);
P_FD=zeros(5001,5);
P_tot=zeros(5001,5);

mm=1;
for n0=[13,14,15,16,17]

% calculation of TA differentiation rate d
syms y
d=solve((1+r/(y-l)+r*y/(g*(y-l)))*n0==ntot,y);
d=double(d);

% calculation of auxiliary quantities 
a=(d-l)/g;
B=d*l/g^2;
rr=r*n0/l;
b=l/g;

h = @(x)(B.*(1 - x)).^(a./2).*besseli(-a, 2.*(B.*(1 - x)).^(1./2));
q1=@(x)a+2*(B*(1-x)).^(1/2).*besseli(a+1,2*(B*(1-x)).^(1/2))./besseli(a,2*(B*(1-x)).^(1/2));
q2=@(x)-a+2*(B*(1-x)).^(1/2).*besseli(-a+1,2*(B*(1-x)).^(1/2))./besseli(-a,2*(B*(1-x)).^(1/2));

% probability-generating function for TA cell population
F_TA=@(x)((d-l)./(d-l.*x)).^rr;
% probability-generating function for FD cell population
F_FD = @(x)(gamma(1 - a).*h(x).*(1 - (a + q2(x))./(a + q1(x)))).^rr;
% probability-generating function for total cell population
F_tot = @(x)x.^(n0).*(gamma(1 - a).*h(x).*(1 - (a +2*b*(1-x)+ q2(x))./(a+2*b*(1-x) + q1(x)))).^rr;


% contour on complex plane
C = [1+1i -1+1i -1-1i 1-1i];

% calculation of probability of TA cell population
parfor k=0:1600
    X1(k+1)=k;
    fun_TA=@(x)F_TA(x)./x.^(k+1);
    % calculation via Cauchy integration
    P_TA(k+1,mm)=real((1/(2*pi*1i))*integral(fun_TA,(1+1i),(1+1i),'Waypoints',C));
end

% calculation of probability of FD and total cell population
parfor k=0:5000
    X2(k+1)=k;
    fun_FD=@(x)F_FD(x)./x.^(k+1);
    fun_tot=@(x)F_tot(x)./x.^(k+1);
    P_FD(k+1,mm)=real((1/(2*pi*1i))*integral(fun_FD,(1+1i),(1+1i),'Waypoints',C));
    P_tot(k+1,mm)=real((1/(2*pi*1i))*integral(fun_tot,(1+1i),(1+1i),'Waypoints',C));
end
mm=mm+1;
end


figure
hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
plot(X1,P_TA(:,1),'color','#9D8678','linewidth',4)
plot(X1,P_TA(:,2),'color','#b09e93','linewidth',4)
plot(X1,P_TA(:,4),'color','#e4bbbf','linewidth',4)
plot(X1,P_TA(:,5),'color','#cdb7a8','linewidth',4)
plot(X1,P_TA(:,3),'color','#a600ff','linewidth',4)
plot(X1,(P_TA(:,1)+P_TA(:,2)+P_TA(:,3)+P_TA(:,4)+P_TA(:,5))/5,':','Color','#00ffa6','linewidth',4)
legend('$N_0=13$','$N_0=14$','$N_0=16$','$N_0=17$','$N_0=15$','averaged distribution','interpreter','latex','fontsize',16)
xlabel('Number of TA cells')
ylabel('Probability mass function')
xlim([0 1600])
ylim([0 2*10^(-3)])

figure
hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
plot(X2,P_FD(:,1),'color','#9D8678','linewidth',4)
plot(X2,P_FD(:,2),'color','#b09e93','linewidth',4)
plot(X2,P_FD(:,4),'color','#e4bbbf','linewidth',4)
plot(X2,P_FD(:,5),'color','#cdb7a8','linewidth',4)
plot(X2,P_FD(:,3),'color','#a600ff','linewidth',4)
plot(X2,(P_FD(:,1)+P_FD(:,2)+P_FD(:,3)+P_FD(:,4)+P_FD(:,5))/5,':','Color','#00ffa6','linewidth',4)
legend('$N_0=13$','$N_0=14$','$N_0=16$','$N_0=17$','$N_0=15$','averaged distribution','interpreter','latex','fontsize',16)
xlabel('Number of FD cells')
ylabel('Probability mass function')
xlim([0 5000])
ylim([0 8*10^(-4)])

figure
hold on
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2;
plot(X2,P_tot(:,1),'color','#9D8678','linewidth',4)
plot(X2,P_tot(:,2),'color','#b09e93','linewidth',4)
plot(X2,P_tot(:,4),'color','#e4bbbf','linewidth',4)
plot(X2,P_tot(:,5),'color','#cdb7a8','linewidth',4)
plot(X2,P_tot(:,3),'color','#a600ff','linewidth',4)
plot(X2,(P_tot(:,1)+P_tot(:,2)+P_tot(:,3)+P_tot(:,4)+P_tot(:,5))/5,':','Color','#00ffa6','linewidth',4)
legend('$N_0=13$','$N_0=14$','$N_0=16$','$N_0=17$','$N_0=15$','averaged distribution','interpreter','latex','fontsize',14)
xlabel('Total number of cell per crypt')
ylabel('Probability mass function')
xlim([0 5000])
ylim([0 6*10^(-4)])
