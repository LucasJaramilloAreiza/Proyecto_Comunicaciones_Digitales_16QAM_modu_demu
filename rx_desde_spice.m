%% rx_desde_spice.m -- canal real + anti-alias + ADC + receptor PSoC sobre la salida de LTspice
% Requiere: salida.txt (LTspice: File > Export data as text, solo V(sout)) y tx_info.mat
clear; clc; close all; rng(3);
load('tx_info.mat');                                   % lI lQ tlead Tsym M Fs fc N
sym = lI + 1j*lQ;  Ts = 1/Fs;
try, D = readmatrix('salida.txt','NumHeaderLines',1);
catch, D = dlmread('salida.txt','',1,0); end
[t,iu] = unique(D(:,1)); s = D(iu,2);
Fh = 16*Fs; R = Fh/Fs;                                 % rejilla fina (4 MS/s)
th = (0:1/Fh:tlead+(N+2)*Tsym).';
sh = interp1(t,s,th,'linear',0);
Ntr = 150;                                             % simbolos de entrenamiento

% --- calibracion con canal ideal: retardo del simbolo y K0 (V por nivel) ---
best = [inf 0 0];
for dl = 0:Tsym/80:0.4*Tsym
    z = zeros(N,1);
    for i = 1:N
        a = round((tlead+(i-1)*Tsym+dl)*Fh)+1; w = a:a+M*R-1;
        z(i) = 2*mean(sh(w).*exp(-2j*pi*fc*th(w)));
    end
    g = (sym'*z)/(sym'*sym); res = norm(z-g*sym)/norm(z);
    if res < best(1), best = [res dl abs(g)]; zcal = z/g; end
end
delay = best(2); K0 = best(3);
fprintf('Retardo del simbolo: %.1f us | K0 = %.3f V/nivel | residuo ideal = %.3f\n',delay*1e6,K0,best(1));

res_ideal = canal_rx(true,  sh, th, delay, K0, lI, lQ, tlead, Tsym, M, Fs, fc, N);
res_real  = canal_rx(false, sh, th, delay, K0, lI, lQ, tlead, Tsym, M, Fs, fc, N);

%% ---------- FIGURA 1: salida del transmisor (SPICE) ----------
figure(1); clf; set(gcf,'Color','w','Name','TX LM1496 (SPICE)');
w = round(tlead*Fh)+(1:4*M*R); subplot(1,3,1); plot((th(w)-th(w(1)))*1e6,sh(w)); grid on;
title('s(t): 4 simbolos'); xlabel('\mus'); ylabel('V');
L = 2^nextpow2(numel(sh)); hw = 0.5-0.5*cos(2*pi*(0:numel(sh)-1).'/(numel(sh)-1));
S = abs(fft(sh.*hw,L)).^2; f = (0:L/2-1)*Fh/L/1e3;
subplot(1,3,2); plot(f,10*log10(S(1:L/2)/max(S)+1e-12)); grid on; xlim([0 10*fc/1e3]); ylim([-80 3]);
title('Espectro (armonicos 3fc, 5fc...)'); xlabel('kHz'); ylabel('dB');
[Lx,Ly] = meshgrid([-3 -1 1 3]);
subplot(1,3,3); plot(real(zcal),imag(zcal),'.'); hold on; plot(Lx(:),Ly(:),'r+'); axis equal; grid on;
title('Constelacion ideal tras LM1496 (compresion/ISI)'); xlim([-5 5]); ylim([-5 5]);

%% ---------- FIGURAS 2 y 3: receptor ----------
graficar(res_ideal,2,'IDEAL'); graficar(res_real,3,'REAL');

%% ================= FUNCIONES LOCALES =================
function Rs = canal_rx(ideal, sh, th, delay, K0, lI, lQ, tlead, Tsym, M, Fs, fc, N)
Ts=1/Fs; sym=lI+1j*lQ; Fh=16*Fs; R=16;
Ntr=150; EsN0dB=25; Pref=10;
al=ones(N,1); if ~ideal, al(201:320)=0.5; al(321:430)=0.3; end
alh=[ones(round(tlead*Fh),1); repelem(al,M*R); ones(numel(th),1)]; alh=alh(1:numel(th));
r=alh.*sh;
if ~ideal
    sig2c=Pref/10^(EsN0dB/10)/2; sn=sqrt(sig2c*K0^2*M/2);
    r=r+sn*sqrt(R)*randn(size(r));                  % misma PSD de ruido a Fh
end
aL=exp(-2*pi*2*fc/Fh);                              % anti-alias 2 polos @ 2fc
r=filter(1-aL,[1 -aL],filter(1-aL,[1 -aL],r));
r=r(1:R:end); tadc=th(1:R:end);                     % ADC a Fs
if ideal, frx=fc; th0=0; else, frx=fc*(1-20e-6); th0=0.6; end   % deriva de 20 ppm
Vfs=4.096; lsb=2*Vfs/4096; pgal=[1 2 4 8 16 24 32 48];
aP=0.05; mu=0.05; zeta=0.707; BnT=0.01;
tt=BnT/(zeta+1/(4*zeta)); d=1+2*zeta*tt+tt^2; Kp=4*zeta*tt/d; Ki=4*tt^2/d;
sl=@(x) min(max(2*floor(x/2)+1,-3),3);
n0=round((tlead+delay)/Ts)+1;
gh=[0 0]; ph=[0 0]; g=0; Pm=Pref; phi=0; integ=0;
z=zeros(N,1); dec=zeros(N,1); G=zeros(N,1); PG=zeros(N,1); PM=zeros(N,1); PH=zeros(N,1); TH=zeros(N,1);
for i=1:N
    gp=gh(1); pp=ph(1);                             % 1 simbolo de latencia en el lazo
    lin=10^(gp/20); pga=max(pgal(pgal<=max(lin,1))); gd=lin/pga;
    idx=n0+(i-1)*M+(0:M-1);
    x=round(min(max(pga*r(idx),-Vfs),Vfs)/lsb)*lsb*gd;
    arg=2*pi*frx*tadc(idx)-th0+pp;
    z(i)=(2/M)*sum(x.*exp(-1j*arg))/K0;
    Pm=(1-aP)*Pm+aP*abs(z(i))^2;
    g=min(max(g+mu*(10*log10(Pref)-10*log10(Pm)),-6),40);
    if i<=Ntr, dd=sym(i); else, dd=sl(real(z(i)))+1j*sl(imag(z(i))); end
    e=imag(z(i)*conj(dd))/abs(dd)^2;
    integ=integ+Ki*e; phi=phi+Kp*e+integ;
    dec(i)=dd; G(i)=g; PG(i)=20*log10(pga); PM(i)=Pm; PH(i)=phi;
    TH(i)=2*pi*(fc-frx)*tadc(idx(1)+M/2)+th0;
    gh=[gh(2) g]; ph=[ph(2) phi];
end
Rs.z=z; Rs.err=(dec~=sym)&((1:N).'>Ntr); Rs.alpha=al; Rs.G=G; Rs.PG=PG; Rs.PM=PM; Rs.PH=PH; Rs.TH=TH;
Rs.ser=movmean(double(Rs.err),40)*100; Rs.Pref=Pref;
fprintf('%s: errores de simbolo = %d / %d  (%.2f %%)\n',ternary(ideal,'IDEAL','REAL '),sum(Rs.err),N-Ntr,mean(Rs.err(Ntr+1:end))*100);
end
function o = ternary(c,a,b), if c, o=a; else, o=b; end, end
function graficar(Rs,fig,nm)
N=numel(Rs.G); k=(0:N-1).'; [Lx,Ly]=meshgrid([-3 -1 1 3]);
figure(fig); clf; set(gcf,'Color','w','Name',['Receptor ' nm]);
subplot(3,2,1); plot(k,20*log10(Rs.alpha),'r',k,Rs.G,'b',k,Rs.PG,'g'); grid on;
legend('\alpha','AGC total','PGA','Location','best'); title(['AGC - ' nm]); xlabel('simbolo'); ylabel('dB');
subplot(3,2,2); plot(k,10*log10(Rs.PM),'b',k,10*log10(Rs.Pref)*ones(N,1),'k--'); grid on; ylim([5 15]); title('P_{med} vs P_{ref}'); xlabel('simbolo'); ylabel('dB');
subplot(3,2,3); plot(k,Rs.TH,'r',k,Rs.PH,'b--'); grid on; legend('\theta real','\theta NCO','Location','best'); title('Costas / NCO'); xlabel('simbolo'); ylabel('rad');
subplot(3,2,4); plot(k,Rs.ser,'m'); grid on; ylim([-1 100]); title('SER (ventana 40) [%]'); xlabel('simbolo');
id=250:320; subplot(3,2,5); plot(real(Rs.z(id)),imag(Rs.z(id)),'.'); hold on; plot(Lx(:),Ly(:),'r+'); axis equal; grid on; xlim([-5 5]); ylim([-5 5]); title('Constelacion (\alpha = 0.5)');
id=450:N; subplot(3,2,6); plot(real(Rs.z(id)),imag(Rs.z(id)),'.'); hold on; plot(Lx(:),Ly(:),'r+'); axis equal; grid on; xlim([-5 5]); ylim([-5 5]); title('Constelacion (canal recuperado)');
end