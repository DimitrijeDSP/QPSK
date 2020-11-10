clear all, close all, clc;

%obrada(ro,trajanje_RRC,prosirenje_NF,fi_0,cf,ck,cs);
%ro : koeficijent zaobljenja
%trajanje_RRC : broj puta trajanja perioda signaliziranja (4,8...)
%prosirenje_NF : koliko je siri u odnosu na original (razlomak)
%fi_0 : fazni pomeraj   (idealno je 0)
%cf : crtaju se rrc i nf filtri samo za cf !=0
%ck : crta se konstelacioni dijagram samo za ck != 0
%cs : crta se sgss samo za cs != 0

% Prikaz filtara i Pe za razne ro i N
obrada(0.25,4,0,0,1,0,0); 
obrada(0.25,8,0,0,1,0,0);   % Ovde se prikazuju i konstelacioni dijagram i sgss
obrada(1,4,0,0,1,0,0);
obrada(1,8,0,0,1,1,1);      % Ovde se prikazuju i konstelacioni dijagram i sgss

% 50% siri NF filtar
obrada(1,8,1/2,0,0,0,0);

% Prikaz konstelacije i Pe za neidealnu sinhronizaciju
obrada(1,8,0,pi/16,0,1,0);