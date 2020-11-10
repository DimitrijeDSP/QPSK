function crtaj_filtar(h)
n = 0:(length(h)-1);
figure,subplot(2,2,1),stem(n,h);
title('Impulsni odziv filtra'),
xlabel('n'),ylabel('h\_[n]'),
xlim([0 (length(h)-1)]);
[H,w] = freqz(h,1,1000);
subplot(2,2,2),plot(w/pi,abs(H));
title('Amplitudska karakteristika filtra'),
xlabel('w/pi'),ylabel('|H\_(e^{jw})|');
subplot(2,2,3),plot(w/pi,unwrap(angle(H)));
title('Fazna karakteristika filtra')
xlabel('w/pi'),ylabel('faza(H\_(e^{jw}))');
subplot(2,2,4),plot(w/pi,20*log10(H),[0 1],[-30 -30],'k--');
title('Pojacanje filtra') 
xlabel('w/pi'),ylabel('20log(H\_(e^{jw}))');