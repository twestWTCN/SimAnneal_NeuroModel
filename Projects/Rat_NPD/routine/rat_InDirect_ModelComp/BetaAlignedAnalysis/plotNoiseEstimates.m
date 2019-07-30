
ptl = linspace(0,100,100)
noiseZ = std(abs(xsim_noise{1}(2,:)));
for i = 1:numel(ptl)

ptileAmp(i) = prctile(abs(xsim_gl{1}(2,:)),ptl(i));

noiseAmp(i) = ptileAmp(i)/noiseZ
end

% subplot(2,1,1)
%  scatter(ptl,ptileAmp)
% xlabel('Percentile of Envelope')
% ylabel('Amp')
% 
% subplot(2,1,2)
 scatter(ptl,10.*log10(noiseAmp))
xlabel('Percentile of Envelope')
ylabel('10 logdB SNR')
ylim([-10 20])

