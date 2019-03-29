function iwave = makeIwave(Iwave_win,Iwave_amp,IwaveStat)
% onset1 = 0.010 + 0.0001*randn;
% onset2 = 0.012 + 0.0001*randn;
onset1 = 0.020; % + 0.0001*randn;
onset2 = 0.022; % + 0.0001*randn;

amp1 = IwaveStat(1,1)*Iwave_amp + IwaveStat(1,2)*randn;
amp2 = IwaveStat(2,1)*Iwave_amp + IwaveStat(2,2)*randn;


[x,i1] = gaussfx(Iwave_win,onset1,0.0005);
i1 = (i1./max(i1)).*amp1;
[x,i2] = gaussfx(Iwave_win,onset2,0.0005);
i2 = (i2./max(i2)).*amp2;

iwave = i1+i2;

% x = x.*1000;
% plot(x,i1,x,i2,x,i1+i2)