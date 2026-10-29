%% SDR 16-QAM A NIVEL DE COMPONENTES CON MC1496/LM1496 (forma de onda, sin toolboxes)
% ESP32 DAC 8b -> RC -> Op-Amp restador -> LM1496 x2 (salidas sumadas en corriente)
% -> canal -> filtro anti-alias -> PSoC
% (PGA discreto + ADC 12b + AGC + NCO/Costas DD). Corre en MATLAB Mobile.
clear; clc; close all;

%% ---------- PARAMETROS DE HARDWARE ----------
P.N=3000; P.Ntr=300; P.M=100;            % simbolos, entrenamiento, muestras/simbolo
P.Fs=500e3; P.fc=50e3; P.Rs=P.Fs/P.M;    % Fs = 10*fc ; Rs = 5 ksym/s
P.Vdd=3.3; P.dacBits=8; P.lsbSym=0.55;   % DAC ESP32: 0.55 V por nivel
P.R=1e3; P.C=10e-9;                      % RC pasivo (fc = 15.9 kHz)
P.Vref=1.65; P.Rin=10e3; P.Rf=10e3;      % restador
% --- LM1496 (par diferencial con degeneracion RE, conmutado por portadora) ---
P.VT=0.026; P.I5=1e-3;                   % V termico, corriente de polarizacion (pin 5)
P.RE=3.3e3; P.RL=1e3; P.Vc=0.5;          % resistencia de ganancia, carga, amplitud portadora [V]
P.re=P.VT/P.I5;                          % re = VT/I5 (cada fuente da I5)
th_=linspace(0,2*pi,20001);
g1=2*mean(tanh(P.Vc*cos(th_)/(2*P.VT)).*cos(th_));   % fundamental de la conmutacion (~4/pi)
P.K0=P.lsbSym*P.Rf/P.Rin*2*P.RL/(P.RE+2*P.re)*g1;   % diferencial = 2*RL/(RE+2re)      % V por nivel (pequena senal)
P.fAA=100e3;                             % anti-alias de 2 polos antes del PGA/ADC
P.Vfs=4.096; P.adcBits=12;               % ADC +-Vfs (ajusta a tu Vref real)
P.pgaSteps=[1 2 4 8 16 24 32 48];        % ganancias PGA (verifica en datasheet)
P.Pref=10; P.aP=0.05; P.muAGC=0.05; P.gmin=-6; P.gmax=40;
P.zeta=0.707; P.BnT=0.01;
P.EsN0dB=25; P.th0=0.6; P.dfHz=1.0;      % canal real: Es/N0, desfase, 1 Hz de deriva

%% ---------- SIMULACION ----------
Ri = simular(P,true);    % canal ideal
Rr = simular(P,false);   % canal real
fprintf('SER ideal: %.2f %% | SER real: %.2f %%\n', ...
    mean(Ri.err(P.Ntr+1:end))*100, mean(Rr.err(P.Ntr+1:end))*100);

%% ---------- FIGURA 1: CADENA ANALOGICA ----------
M=P.M; Fs=P.Fs; w=1:8*M; tu=(w-1)/Fs*1e6; ws=1:4*M;
figure(1); clf; set(gcf,'Color','w','Name','Cadena Tx');
subplot(3,2,1); stairs(tu,Ri.vdacI(w),'r'); hold on; plot(tu,Ri.vrcI(w),'b','LineWidth',1.2);
grid on; legend('DAC 8b','tras RC'); title('ESP32 DAC (I) y filtro RC'); xlabel('\mus'); ylabel('V');
subplot(3,2,2); plot(tu,Ri.vI(w),'b',tu,Ri.vQ(w),'r'); grid on;
legend('I bipolar','Q bipolar'); title('Salida op-amp restador'); xlabel('\mus'); ylabel('V');
subplot(3,2,3); plot(tu(ws),Ri.s(ws),'k'); grid on;
title('s(t) LM1496 (I+Q en corriente)'); xlabel('\mus'); ylabel('V');
c=(M/2):M:P.N*M; subplot(3,2,4);
plot(Ri.vI(c),Ri.vQ(c),'b.','MarkerSize',3); axis equal; grid on;
title('Constelacion X-Y (osciloscopio)'); xlabel('I [V]'); ylabel('Q [V]');
L=2^nextpow2(numel(Ri.s)); win=0.5-0.5*cos(2*pi*(0:numel(Ri.s)-1).'/(numel(Ri.s)-1));
S=abs(fft(Ri.s.*win,L)).^2; f=(0:L/2-1)*Fs/L/1e3;
subplot(3,2,5); plot(f,10*log10(S(1:L/2)/max(S)+1e-12)); grid on; xlim([0 150]); ylim([-80 3]);
title('Espectro de s(t)'); xlabel('kHz'); ylabel('dB');
w2=740*M+(1:20*M); subplot(3,2,6);
plot((w2-w2(1))/Fs*1e3,Rr.r(w2),'r',(w2-w2(1))/Fs*1e3,Ri.r(w2),'b'); grid on;
legend('canal real','ideal'); title('r(t): escalon de atenuacion'); xlabel('ms'); ylabel('V');

%% ---------- FIGURAS 2 y 3: RECEPTOR PSoC ----------
graficar(Rr,P,'RECEPTOR - CANAL REAL',2);
graficar(Ri,P,'RECEPTOR - CANAL IDEAL',3);

%% ================= FUNCIONES LOCALES =================
function R = simular(P,ideal)
    rng(1); N=P.N; M=P.M; Fs=P.Fs; fc=P.fc;
    b=randi([0 1],N,4);
    lv=@(b1,b0) 2*(2*b1+mod(b1+b0,2))-3;
    lI=lv(b(:,1),b(:,2)); lQ=lv(b(:,3),b(:,4)); sym=lI+1j*lQ;
    % ESP32 DAC (cuantizado)
    q=@(v) round(min(max(v,0),P.Vdd)/P.Vdd*(2^P.dacBits-1))/(2^P.dacBits-1)*P.Vdd;
    R.vdacI=q(repelem(P.Vref+P.lsbSym*lI,M,1));
    R.vdacQ=q(repelem(P.Vref+P.lsbSym*lQ,M,1));
    % Filtro RC
    a=exp(-1/(Fs*P.R*P.C));
    R.vrcI=filter(1-a,[1 -a],R.vdacI,P.Vref*a);
    R.vrcQ=filter(1-a,[1 -a],R.vdacQ,P.Vref*a);
    % Op-amp restador (tolerancia del divisor y offset en el caso real)
    if ideal, tol=0; Vos=0; else, tol=0.01; Vos=5e-3; end
    Vr=P.Vref*(1+tol);
    R.vI=P.Rf/P.Rin*(R.vrcI-Vr)+Vos; R.vQ=P.Rf/P.Rin*(R.vrcQ-Vr)+Vos;
    % LM1496 x2: i = I5*tanh(Vs/(I5*(RE+2re))) * tanh(Vc/(2VT)); salidas en paralelo -> suma de corrientes
    t=(0:N*M-1).'/Fs;
    dI=@(v) 2*P.I5*tanh(v/(P.I5*(P.RE+2*P.re)));   % corriente diferencial, satura en 2*I5
    sw=@(c) tanh(c/(2*P.VT));
    R.s=P.RL*( dI(R.vI).*sw(P.Vc*cos(2*pi*fc*t)) + dI(R.vQ).*sw(-P.Vc*sin(2*pi*fc*t)) );
    % Canal: atenuador + AWGN
    al=ones(N,1);
    if ~ideal, al(751:1750)=0.5; al(1751:2500)=0.3; end
    if ideal, sn=0; else
        sig2c=P.Pref/10^(P.EsN0dB/10)/2; sn=sqrt(sig2c*P.K0^2*M/2);
    end
    R.r=repelem(al,M,1).*R.s + sn*randn(N*M,1);
    aL=exp(-2*pi*P.fAA/Fs);                                  % anti-alias (quita armonicos 9fc, 11fc...)
    R.r=filter(1-aL,[1 -aL],filter(1-aL,[1 -aL],R.r));
    R.alpha=al;
    if ideal, frx=fc; th0=0; else, frx=fc-P.dfHz; th0=P.th0; end
    % Receptor PSoC
    tt=P.BnT/(P.zeta+1/(4*P.zeta)); d=1+2*P.zeta*tt+tt^2;
    Kp=4*P.zeta*tt/d; Ki=4*tt^2/d;
    sl=@(x) min(max(2*floor(x/2)+1,-3),3);
    lsb=2*P.Vfs/2^P.adcBits;
    g=0; Pm=P.Pref; phi=0; integ=0;
    z=zeros(N,1); dec=zeros(N,1); R.gdB=zeros(N,1); R.pgadB=zeros(N,1);
    R.Pm=zeros(N,1); R.phi=zeros(N,1);
    for i=1:N
        idx=(i-1)*M+(1:M);
        lin=10^(g/20); pga=max(P.pgaSteps(P.pgaSteps<=max(lin,1))); gd=lin/pga;
        x=round(min(max(pga*R.r(idx),-P.Vfs),P.Vfs)/lsb)*lsb*gd;   % PGA->ADC->ganancia digital
        arg=2*pi*frx*t(idx)-th0+phi;                               % NCO (LUT)
        z(i)=(2/M)*sum(x.*exp(-1j*arg))/P.K0;                      % I/Q (integrate&dump)
        Pm=(1-P.aP)*Pm+P.aP*abs(z(i))^2;
        g=min(max(g+P.muAGC*(10*log10(P.Pref)-10*log10(Pm)),P.gmin),P.gmax);
        if i<=P.Ntr, dd=sym(i); else, dd=sl(real(z(i)))+1j*sl(imag(z(i))); end
        e=imag(z(i)*conj(dd))/abs(dd)^2;
        integ=integ+Ki*e; phi=phi+Kp*e+integ;
        dec(i)=dd; R.gdB(i)=g; R.pgadB(i)=20*log10(pga); R.Pm(i)=Pm; R.phi(i)=phi;
    end
    R.z=z; R.err=(dec~=sym)&((1:N).'>P.Ntr);
    R.ser=movmean(double(R.err),200)*100;
    tc=((0:N-1).'*M+M/2)/Fs;
    R.thTrue=2*pi*(fc-frx)*tc+th0;
end

function graficar(R,P,ttl,fig)
    N=P.N; k=(0:N-1).';
    figure(fig); clf; set(gcf,'Color','w','Name',ttl);
    subplot(3,2,1); plot(k,20*log10(R.alpha),'r',k,R.gdB,'b',k,R.pgadB,'g'); grid on;
    legend('\alpha canal','AGC total','PGA (analogico)','Location','best');
    title([ttl ': AGC']); xlabel('simbolo'); ylabel('dB');
    subplot(3,2,2); plot(k,10*log10(R.Pm),'b',k,10*log10(P.Pref)*ones(N,1),'k--'); grid on;
    ylim([5 15]); title('Potencia vs P_{ref}'); xlabel('simbolo'); ylabel('dB');
    subplot(3,2,3); plot(k,R.thTrue,'r',k,R.phi,'b--'); grid on;
    legend('\theta real','\theta NCO','Location','best'); title('Costas / NCO'); xlabel('simbolo'); ylabel('rad');
    subplot(3,2,4); plot(k,R.ser,'m'); grid on; ylim([-1 15]);
    title('SER (ventana 200)'); xlabel('simbolo'); ylabel('%');
    id=1000:1700; subplot(3,2,5);
    plot(real(R.z(id)),imag(R.z(id)),'.','MarkerSize',4); hold on;
    [Lx,Ly]=meshgrid([-3 -1 1 3]); plot(Lx(:),Ly(:),'r+','MarkerSize',8);
    axis equal; grid on; xlim([-5 5]); ylim([-5 5]); title('Constelacion (alpha=0.5)');
    id=2600:3000; subplot(3,2,6);
    plot(real(R.z(id)),imag(R.z(id)),'.','MarkerSize',4); hold on; plot(Lx(:),Ly(:),'r+','MarkerSize',8);
    axis equal; grid on; xlim([-5 5]); ylim([-5 5]); title('Constelacion (canal recuperado)');
end
