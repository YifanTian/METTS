
coeificient = load('coefficient.txt');

E_1 = load('tebd_energybond.txt');
E_2 = load('tebd_energybond_sbs.txt');
E_3 = zeros(1,59);
for i = 1:59
    E_3(i) = 1/4-log(2);
end
txi = load('coefficient.txt'); 
lattice = [1:1:59];
plot(lattice,[E_1(:),E_2(:),E_3(:)],'LineWidth',1.8)
title('TEBD Heisenberg spin-half, energy bond, T = 10');
legend('open','smooth','L = infinity,J=1.0,L=60')
xlabel('i')
ylabel('S_iS_{i+1}')
set(gca,'fontsize',16)


x = [0:1:40];
c1y = zeros(1,40);
c2y = zeros(1,40);
tx = [1:1:60];
for i = 1:length(x)
    c1y(i) = c1(x(i),40);
    c2y(i) = c2(x(i),20);
end

% plot(x,c1y(:),'LineWidth',1.8)
% hold on
% plot(x,c2y(:),'o','LineWidth',1.8)
% legend('y(m/M)','Borel')
% xlabel('m')
% ylabel('c_m')
% set(gca,'fontsize',16)

M = 20;
L = length(tx);
ti = zeros(1,L);
for i = 1:L
    if (i >= 1) && (i <= M)
        ti(i) = c1(M-i,M);
    elseif (i>M) && (i <= L-M)
        %fprintf('here1 %d\n',i);
        ti(i) = 1;
    elseif (i > L-M) && (i<=L)
        %fprintf('here2 %d %d\n ',i,L);
        ti(i) = c1(i-L+M,M);
    end
end


plot(tx,ti,'LineWidth',1.8)
title('smooth boundary condition,J_i/J');
xlabel('i')
ylabel('J_i/J')
set(gca,'fontsize',16)

function f = c1(m,M)
    x = m/M;
    f = 0.5*(1-tanh((x-0.5)/(x*(1-x))));
end

function f = c2(m,M)
    curr = 0;
    for n = m+1:50
        curr = curr + M^n/factorial(n);
    end
    f = curr*exp(-M);
end



