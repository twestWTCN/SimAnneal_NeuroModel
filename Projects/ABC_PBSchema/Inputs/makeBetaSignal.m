function [t,tx,px,ax] = makeBetaSignal(t_in)
fsamp = 1./(t_in(5)-t_in(4));
% Create background beta oscillation
lorentzian = @(x,kx,gamma)((1/(pi*gamma))*((gamma^2)./((x-kx).^2 + (gamma^2))));
% x = 0:fsamp/N:fsamp/2;
x = 0:0.001:(fsamp); % support
kx = 20; % centre
gamma = 0.5; %FWHM
fx = lorentzian(x,kx,gamma) .*exp(1i.*rand(size(x)));
% figure
% subplot(2,1,1)
% plot(x,abs(fx));
% xlim([0 50]);
% xlabel('Hz'); ylabel('Amplitude')

tx = ifft(fx);
tx = real(tx(fsamp+1:end-fsamp-1));

t = linspace(0,size(tx,2)./(fsamp),size(tx,2));
trncval = find(t<t_in(end),1,'last');
t = t(1:trncval);
tx = tx(1:trncval);
tx = (tx-mean(tx))./std(tx);

px = wrapToPi(angle(hilbert(tx)));
ax = abs(hilbert(tx));

% subplot(2,1,2)
% plot(t,tx);
% xlabel('Time'); ylabel('Amplitude')
% xlim([5 6])
% a= 1;