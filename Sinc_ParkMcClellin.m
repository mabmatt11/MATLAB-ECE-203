L = 100;
fm = 1/8;

k = 0:2*L;
h = sin(2*pi*fm*(k-L))./(pi*(k-L));
h(L+1) = 2*fm;

figure(1); 
clf;
stem(k,h,'linewidth',3);
xlabel('k');
ylabel('h[k]');
title(['truncated sinc impulse response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
set(gca,'fontsize',18)

H = fft([h, zeros(1,800)]);
f = linspace(-1/2,1/2,length(H));
H_ideal = 1.*(abs(f)<=fm);
figure(2);plot(f,[abs(fftshift(H));H_ideal],'linewidth',3);
xlabel('$\hat f$','Interpreter','latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title(['truncated sinc frequency response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
legend('truncated sinc','ideal')
set(gca,'fontsize',18)

figure(3);semilogy(f,[abs(fftshift(H))],'linewidth',3);
xlabel('$\hat f$','Interpreter','latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title(['truncated sinc frequency response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
set(gca,'fontsize',18)

figure(4);semilogy(f,[abs(fftshift(H));H_ideal+eps],'linewidth',3);
xlabel('$\hat f$','Interpreter','latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title(['truncated sinc frequency response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
legend('truncated sinc',['ideal + ' num2str(eps)])
set(gca,'fontsize',18)

%%
delta = 1/80;
h_pm = firpm(2*L,[0 fm fm+delta .5]*2,[1 1 0 0]);

figure(11); 
clf;
stem(k,[h;h_pm]','linewidth',3);
xlabel('k');
ylabel('h[k]');
legend('trucated sinc','parks-mcclennan');
title(['truncated sinc impulse response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
set(gca,'fontsize',18)

H_pm = fft([h_pm, zeros(1,800)]);
f = linspace(-1/2,1/2,length(H_pm));
figure(12);plot(f,[abs(fftshift(H));abs(fftshift(H_pm));H_ideal],'linewidth',3);
xlabel('$\hat f$','Interpreter','latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title(['truncated sinc frequency response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
legend('trucated sinc','parks-mcclennan','ideal');
set(gca,'fontsize',18)

figure(13);semilogy(f,[abs(fftshift(H));abs(fftshift(H_pm))],'linewidth',3);
xlabel('$\hat f$','Interpreter','latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title(['truncated sinc frequency response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
legend('trucated sinc','parks-mcclennan');
set(gca,'fontsize',18)

figure(14);semilogy(f,[abs(fftshift(H));abs(fftshift(H_pm));H_ideal+eps],'linewidth',3);
xlabel('$\hat f$','Interpreter','latex');
ylabel('$|H(\hat f)|$','Interpreter','latex');
title(['truncated sinc frequency response, $f_m =$' num2str(fm) ', $L =$' num2str(L)],'Interpreter','latex');
legend('trucated sinc','parks-mcclennan','ideal');
set(gca,'fontsize',18)