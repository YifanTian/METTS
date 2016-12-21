
ensemble = load('ensemble.txt');
plot(ensemble(:,1),ensemble(:,2),'-','LineWidth',1.8)

E_thermal = load('thermal_step.txt');
plot(E_thermal(:,1),E_thermal(:,2),'-','LineWidth',1.8)


E2_1 = load('Metts_FE.txt'); 
E1_1 = load('FiniteTE.txt'); 
%plot(E2_1(:,1),E2_1(:,2),'o','LineWidth',1.8)
%errorbar(E2_1(:,1),E2_1(:,2),E2_1(:,3))
errorbar(E2_1(:,1),E2_1(:,2),E2_1(:,3),'-s','MarkerSize',2,...
    'MarkerEdgeColor','black','MarkerFaceColor','black');
grid on
hold on
plot(E1_1(:,1),E1_1(:,2),'-','LineWidth',1.0,'MarkerEdgeColor','y',...
     'MarkerFaceColor',[0.5,0.5,0.5])
legend('METTS','exact finite')
%legend('exact finite','exact ground','ancilla')
title('Finite temperature, Heisenberg spin-half');
xlabel('Temperature(beta)')
ylabel('Energy')
set(gca,'fontsize',16)


%==========================================
E1_1 = load('FiniteTE.txt'); 
E_finite = ones(1000,1)*(E1_1(40,2));
E_exact = ones(1000,1)*(-4.2580352);
Et_2 = load('thermal_step_2.txt');
%plot(Et_2(:,1),Et_2(:,2),'o','LineWidth',1.0)
%plot(Et_2(:,1),[Et_2(:,2),Et_2(:,3),E_exact],'o','LineWidth',1.0)
plot(Et_2(:,1),Et_2(:,2),'o','LineWidth',1.0)
hold on
plot(Et_2(:,1),Et_2(:,3),'-','LineWidth',3.0)
hold on
plot(Et_2(:,1),E_finite,'-','LineWidth',3.0)
hold on
plot(Et_2(:,1),E_exact,'-','LineWidth',3.0)
legend('energy sample','energy average','exact  finite','exact ground')
%legend('energy sample','energy average')
xlabel('thermal step')
ylabel('Energy')
title('Finite temperature beta = 8, Heisenberg spin-half');
set(gca,'fontsize',16)
%plot(Et_2(:,1),Et_2(:,4),'o','LineWidth',1.8)


Et_8 = load('thermal_step_8.txt');
%plot(Et_8(:,1),Et_8(:,2),'o','LineWidth',1.8)

plot(Et_8(:,1),[Et_8(:,2),Et_8(:,3)],'o','LineWidth',1.8)
plot(Et_8(:,1),Et_8(:,4),'o','LineWidth',1.8)

%==========================================
E2_1 = load('Metts_FE.txt'); 
E1_1 = load('FiniteTE_sz0.txt'); 
E1_2 = load('FiniteTE.txt'); 
plot(E1_1(:,1),E1_2(:,2),'g','LineWidth',1.0)
hold on
plot(E1_1(:,1),E1_1(:,2),'r','LineWidth',1.0)
%hold on
%plot(E1_1(:,1),[E1_1(:,2),E1_2(:,2)],'-','LineWidth',1.0,'MarkerEdgeColor','b',...
%      'MarkerFaceColor',[0.5,0.5,0.5])
hold on
errorbar(E2_1(:,1),E2_1(:,2),E2_1(:,3),'-s','MarkerSize',2,...
    'MarkerEdgeColor','black','MarkerFaceColor','green','LineWidth',1.0);
legend('exact finite','exact finite sz0','METTS')
%legend('exact finite','exact ground','ancilla')
title('Finite temperature, Heisenberg spin-half');
xlabel('Temperature(beta)')
ylabel('Energy')
set(gca,'fontsize',16)
%==========================================

%E20 = load('myfile_1_20.txt'); 
%figure(1); plot(E2_1(:,1),[E1(:,2),E2_1(:,2),E3(:,1)],'-o','LineWidth',1.5);
figure(1); plotfit(E1(:,1),[E1(:,2),E3(:,1)],'-o','LineWidth',1.5);
hold on
figure(1); plot(E2_1(:,1),[E2_1(:,2),E3(:,1)],'-o','LineWidth',1.5);
%legend('termal density matrix','ancilla','exact')
%ylim([-45,-35])
%ylim([-0.5,2])
grid on
%legend('energy')
title('Finite temperature');
xlabel('temperature')
ylabel('energy')


% Sn = load('myfile_Sn.txt');
% %C = load('myfile_C.txt');
% SS = load('myfile_SS.txt');
% figure(2); plot(Sn(1:99,1),Sn(1:99,2:3),'-','LineWidth',1.5);
% legend('Sz','Sx')
% hold on
% figure(2); plot(Sn(1:99,1),SS(:,2),'o','LineWidth',1.5);

%legend('SS')