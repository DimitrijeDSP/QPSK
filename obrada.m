function obrada(ro,trajanje_RRC,prosirenje_NF,fi_0,cf,ck,cs)

pkg load communications;
pkg load signal;

% Protok simbola i binarni protok
M = 2^2;                    % QPSK
V_simbola = 2000;           %sim/s
V_b = log2(M)*V_simbola;    % Ekvivalentni binarni protok

% Binarna pilot sekvenca za frejm
bin_pilot_duz = 60;   % broj pilot bita
bin_pilot = (((-1).^(0:bin_pilot_duz-1)) + 1) / 2;  % Naizmenicne 1 i 0

% Mapiranje pilot bita u kompleksne QPSK simbole
sim_pilot_duz = bin_pilot_duz/log2(M);
sim_pilot = zeros(1,sim_pilot_duz);    % 30 pilot simbola

for br = 1:2:(bin_pilot_duz-1)
  temp = bin_pilot(br:br+1);
  if (temp(1) == 1) && (temp(2) == 1)
    faza = pi/4;
  elseif (temp(1) == 1) && (temp(2) == 0)
    faza = 3*pi/4;
  elseif (temp(1) == 0) && (temp(2) == 0)
    faza = 5*pi/4;
  elseif (temp(1) == 0) && (temp(2) == 1)
    faza = 7*pi/4;
  endif
  sim_pilot(1,(br+1)/2) = e^(i*faza);
endfor

% Generisanje binarne (informacione) sekvence za frejm
P_0 = 1/2;
P_1 = 1/2;
N_bita = 56000;
x = randsrc(N_bita,1,[1 0; P_1 P_0]);

% Mapiranje bita u kompleksne simbole za QPSK
sim_info_duz = N_bita/log2(M);
sim_info = zeros(1,sim_info_duz);

for br = 1:2:(N_bita-1)
  temp = x(br:br+1);
  if (temp(1) == 1) && (temp(2) == 1)
    faza = pi/4;
  elseif (temp(1) == 1) && (temp(2) == 0)
    faza = 3*pi/4;
  elseif (temp(1) == 0) && (temp(2) == 0)
    faza = 5*pi/4;
  elseif (temp(1) == 0) && (temp(2) == 1)
    faza = 7*pi/4;
  endif
  sim_info(1,(br+1)/2) = e^(i*faza);
endfor

% Formiranje frejma sa 170 simbola
info_sim_po_frejmu = 140;                     
f_duz = sim_pilot_duz + info_sim_po_frejmu;   % =170

N_frejmova = sim_info_duz/info_sim_po_frejmu;  % =200

frejm = zeros(N_frejmova, f_duz);     % 200 frejmova po vrstama

for br = 1:N_frejmova
  frejm(br,1:sim_pilot_duz) = sim_pilot;
  frejm(br,(sim_pilot_duz+1):f_duz) = sim_info(1,(br-1) * info_sim_po_frejmu + (1:info_sim_po_frejmu));
endfor

sim_info_rx = zeros(1,sim_info_duz);

Pe_b = zeros(8,1);
Pe_b_teorija = zeros(8,1);

M1up = 32;
M2up = 32;

duz_imp = trajanje_RRC*M1up; 
N_RRC = duz_imp - 1;    % Red filtra
if mod(N_RRC,2) != 0
   N_RRC = N_RRC + 1;
endif
kasnjenje_RRC = N_RRC/2;

% Za duzinu impulsnog odziva 8 se tacno dobija minimalno slabljenje od 30dB
h_RRC = srrcf(N_RRC,M1up,ro);
if cf != 0
  crtaj_filtar(h_RRC);       % Moja funkcija
endif
      
w_us = pi/M2up;     % granicna ucestanost NF filtra
w_us = w_us*(1+prosirenje_NF);
N_NF = ceil(8*pi/w_us);   % PO TEORIJI
if mod(N_NF,2) != 0
   N_NF = N_NF + 1;
endif
kasnjenje_NF = N_NF/2;    
h_NF = fir1(N_NF,w_us/pi,hann(N_NF+1));     % Hanov prozor ima minimalno slabljenje 44 dB u nepropusno opsegu
if cf != 0
  crtaj_filtar(h_NF); 
endif
    
tacke = [e^(i*pi/4) e^(i*3*pi/4) e^(i*5*pi/4) e^(i*7*pi/4)];
N_fft = 4096;   % Duzina podniza za estimaciju spektra

% POSTO SE DUGO IZVRSAVA, za konkretnu vrednost Eb_pN_dB:
% STAVITI for Eb_pN_dB = 7 UMESTO OVOG for ISPOD (ili neku drugu vrednost)
for Eb_pN_dB = [0:7 15] % dB
  
  for k = 1:N_frejmova
    
    % Razdvajanje na granu u fazi i kvadraturi
    u_I_pocetni = real(frejm(k,:));
    u_Q_pocetni = imag(frejm(k,:));
    
    f_sampling = V_simbola;
    
    % Upsampling - faktor M1up
    
    
    u_I_1up = M1up*upsample(u_I_pocetni,M1up);  % Mnozenje sa M1up zbog upsample
    u_Q_1up = M1up*upsample(u_Q_pocetni,M1up);
    
    f_sampling = f_sampling * M1up;
    
    % Uoblicavanje RRC filtrom
    u_I_1up_temp = [u_I_1up zeros(1,kasnjenje_RRC)];
    u_I_1up_RRC_temp = filter(h_RRC,1,u_I_1up_temp);
    u_I_1up_RRC = u_I_1up_RRC_temp((1+kasnjenje_RRC):end);
    
    u_Q_1up_temp = [u_Q_1up zeros(1,kasnjenje_RRC)];
    u_Q_1up_RRC_temp = filter(h_RRC,1,u_Q_1up_temp);
    u_Q_1up_RRC = u_Q_1up_RRC_temp((1+kasnjenje_RRC):end);
    
    % Upsampling - faktor M2up
    u_I_2up = M2up*upsample(u_I_1up_RRC,M2up); % Mnozenje sa M2up zbog upsample
    u_Q_2up = M2up*upsample(u_Q_1up_RRC,M2up);
    
    f_sampling = f_sampling * M2up;   % Sistemska ucestanost odabiranja = V_simbola*M1up*M2up
    
    % Filtriranje ekvivalentnim NF filtrom
    u_I_2up_temp = [u_I_2up zeros(1,kasnjenje_NF)];   % Prosirivanje za duzinu kasnjenja (visak)
    u_I_temp = filter(h_NF,1,u_I_2up_temp);  % Unosi se kasnjenje
    u_I = u_I_temp((1+kasnjenje_NF):end);   % Skidanje viska sa pocetka
    
    u_Q_2up_temp = [u_Q_2up zeros(1,kasnjenje_NF)];
    u_Q_temp = filter(h_NF,1,u_Q_2up_temp);
    u_Q = u_Q_temp((1+kasnjenje_NF):end);
    
    % Mnozenje sa nosiocem u fazi i kvadraturi
    f_0 = f_sampling / 4; % Ucestanost nosioca
    t=0:(1/f_sampling):(length(u_I)-1)/f_sampling;
    %t = [0:(length(u_I)-1)]/f_sampling;
    u_I_mod = u_I.*cos(2*pi*f_0*t);
    u_Q_mod = u_Q.*sin(2*pi*f_0*t);
    
    % Sabiranje, izlaz predajnika, ulaz prijemnika
    u_QPSK = u_I_mod + u_Q_mod;
    
    % Dodavanje AWGN
    BitsBySymbol = log2(M); % =2
    SampleBySymbol = f_sampling / V_simbola;    % =1024 = M1up*M2up
    ShapingFactor = ro;
    SNRdB = Eb_pN_dB + 10*log10(BitsBySymbol) - 10*log10(SampleBySymbol/(1+ShapingFactor));
    u_r = awgn(u_QPSK, SNRdB, 'measured');
    %u_r = u_QPSK; %OVO JE BEZ SUMA
    
    % Mnozenje sa nosiocem u fazi i kvadraturi
    t=0:1/f_sampling:(length(u_r)-1)/f_sampling;
    u_I_t_demod = 2*u_r.*cos(2*pi*f_0*t + fi_0);  % Mnozenje sa 2 zbog proizvoda kosinusa
    u_Q_t_demod = 2*u_r.*sin(2*pi*f_0*t + fi_0);
    
    % Filtriranje ekvivalentnim NF filtrom
    u_I_t_demod_temp = [u_I_t_demod zeros(1,kasnjenje_NF)];
    u_I_t_NF_temp = filter(h_NF,1,u_I_t_demod_temp);
    u_I_t_NF = u_I_t_NF_temp((1+kasnjenje_NF):end);
    
    u_Q_t_demod_temp = [u_Q_t_demod zeros(1,kasnjenje_NF)];
    u_Q_t_NF_temp = filter(h_NF,1,u_Q_t_demod_temp);
    u_Q_t_NF = u_Q_t_NF_temp((1+kasnjenje_NF):end);
    
    % Downsampling - faktor M2up
    u_I_1dn = downsample(u_I_t_NF,M2up);
    u_Q_1dn = downsample(u_Q_t_NF,M2up);
    %f_sampling = f_sampling / M2up; 
    
    % Optimalno filtriranje RRC filtrom
    u_I_1dn_temp = [u_I_1dn zeros(1,kasnjenje_RRC)];
    u_I_t_RRC_temp = filter(h_RRC,1,u_I_1dn_temp);
    u_I_t_RRC = u_I_t_RRC_temp((1+kasnjenje_RRC):end);
    
    u_Q_1dn_temp = [u_Q_1dn zeros(1,kasnjenje_RRC)];
    u_Q_t_RRC_temp = filter(h_RRC,1,u_Q_1dn_temp);
    u_Q_t_RRC = u_Q_t_RRC_temp((1+kasnjenje_RRC):end);
    
    % Downsampling - faktor M1up
    u_I_2dn = downsample(u_I_t_RRC,M2up);
    u_Q_2dn = downsample(u_Q_t_RRC,M2up);
    %f_sampling = f_sampling / M1up;
    
    % Formiranje kompleksnih simbola
    u_t = u_I_2dn + i*u_Q_2dn;
    
    if k <= 3
      prva3f((k-1)*length(u_t) + (1:length(u_t))) = u_t;    % prva 3 frejma za dijagram
    endif
      
    % Odlucivanje : kombinacija znaka realnog i imaginarnog dela odgovara kvadrantu
    frejm_prijem = (sign(real(u_t)) + i*sign(imag(u_t)))/sqrt(2);
    
    % Izdvajanje informacionih simbola 
    sim_info_rx(1,(1:info_sim_po_frejmu)+(k-1)*info_sim_po_frejmu) = frejm_prijem(1,(sim_pilot_duz+1):f_duz);
  
  endfor  % k

  % Demapiranje simbola u bite
  y = zeros(N_bita,1);
  dva_bita = zeros(1,2);

  for k = 1:sim_info_duz
    re = real(sim_info_rx(k));
    im = imag(sim_info_rx(k));
    if (re > 0) && (im > 0)
      dva_bita = [1 1];
    elseif (re < 0) && (im > 0)
      dva_bita = [1 0];
    elseif (re < 0) && (im < 0)
      dva_bita = [0 0];
    elseif (re > 0) && (im < 0)
      dva_bita = [0 1];  
    endif
    y(2*(k-1)+(1:2),1) = dva_bita;
  endfor

  if Eb_pN_dB != 15
    Eb_pN_dB
    Pe_b(Eb_pN_dB+1) = mean(x~=y)
    Eb_pN = 10^(Eb_pN_dB/10);
    Pe_b_teorija(Eb_pN_dB+1) = 0.5*erfc(sqrt(Eb_pN))    % QPSK
  endif
  
  % Konstelacioni dijagram
  if (Eb_pN_dB == 15) && (ck != 0)   
    figure,plot(real(prva3f),imag(prva3f),'bo');
    grid on
    hold on
    plot(real(tacke),imag(tacke),'rx');
    title('Konstelacioni dijagram QPSK'),
    legend('Dobijene vrednosti','Idealne vrednosti'),
    xlabel('komponenta u fazi'),ylabel('komponenta u kvadraturi');
   endif
  
endfor    % Eb_pN_dB

% Estimacija spektra
if cs != 0
  procena_sgss(u_r,N_fft,f_sampling);       % (moja funkcija)
endif

N_sum = [0:7];
% Prikaz verovatnoce greske za Eb_pN_dB = [0 7]
figure,plot(N_sum,[log10(Pe_b) log10(Pe_b_teorija)]),
xlabel('Eb/pN [dB]'),ylabel('log_{10}Peb'),
title('Estimacija verovatnoce greske'),
legend('Procenjena','Teorijska');

