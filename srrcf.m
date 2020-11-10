function h=srrcf(N,M,r)
% N - red filtra
% M - faktor promene u?estanosti
% r - roll-off faktor
n=0:N;
n1=-N/2:N/2;
h=(4*r*(n-N/2).*cos(pi*(n-N/2)*(1+r)/M)+M*sin(pi*(n-N/2)*(1-r)/M))./((1-(4*r*(n-N/2)/M).^2)*pi.*(n-N/2)*M);
h((n-N/2)==M/4/r)=-r/M*(2/pi*cos(pi/4/r*(1+r))-cos(pi/4/r*(1-r)));
h((n-N/2)==-M/4/r)=-r/M*(2/pi*cos(pi/4/r*(1+r))-cos(pi/4/r*(1-r)));
h((n-N/2)==0)=1/M+r/M*(4/pi-1);