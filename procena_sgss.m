function procena_sgss(u_r,N_fft,f_sampling)
N_podn = floor(length(u_r)/N_fft);
podn_u_r = zeros(1,N_fft);
spektar_podn = zeros(N_podn,N_fft);
sgss_podn = zeros(N_podn,N_fft);


for i_podn = 1:N_podn
podn_u_r(1,:) = u_r(1,(i_podn-1)*N_fft+(1:N_fft));
spektar_podn(i_podn,:) = fft(podn_u_r(1,:));
sgss_podn(i_podn,:) = abs(spektar_podn(i_podn,:)).^2 / (f_sampling*N_fft);
endfor

ukupna_sgss = sum(sgss_podn)/N_fft;   
f_dft = [0:N_fft-1]/N_fft*f_sampling;
figure,plot(f_dft(1:(N_fft/2)),10*log10(ukupna_sgss(1:N_fft/2))),
title('Ukupna procena SGSS modulisanog signala na prijemu'),
xlabel('f [Hz]'),ylabel('SGSS [dB/Hz]'),
xlim([0 f_dft(N_fft/2)]);