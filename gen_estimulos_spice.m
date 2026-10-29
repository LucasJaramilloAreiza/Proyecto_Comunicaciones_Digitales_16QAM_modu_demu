%% gen_estimulos_spice.m -- genera los archivos que lee tx_lm1496.cir (LTspice)
% Salida: vdac_I.pwl, vdac_Q.pwl, params.inc, tx_info.mat
clear; clc; rng(1);
fc = 25e3;  M = 40;  N = 500;                 % portadora, muestras/simbolo, simbolos
Fs = 10*fc; Ts = 1/Fs; Tsym = M*Ts;           % Fs = 250 kS/s, Rs = 6.25 ksym/s
tlead = 8/fc;                                 % 8 ciclos de portadora antes del primer simbolo

% bits -> niveles (Gray): 00->-3 01->-1 11->+1 10->+3
b  = randi([0 1], N, 4);
lv = @(b1,b0) 2*(2*b1 + mod(b1+b0,2)) - 3;
lI = lv(b(:,1),b(:,2));  lQ = lv(b(:,3),b(:,4));

% DAC ESP32: 0..3.3 V, 8 bits, 0.55 V por nivel centrado en 1.65 V
q8 = @(v) round(min(max(v,0),3.3)/3.3*255)/255*3.3;
escribir_pwl('vdac_I.pwl', q8(1.65+0.55*lI), tlead, Tsym);
escribir_pwl('vdac_Q.pwl', q8(1.65+0.55*lQ), tlead, Tsym);

Tstop = tlead + N*Tsym + 2*Tsym;
fid = fopen('params.inc','w');
fprintf(fid,'.param fc=%g\n.param Tc={1/fc}\n.param Vc=0.5\n.param Vcm=3\n',fc);
fprintf(fid,'.param RE=3.3k\n.param RL=1k\n.param R5=6.8k\n.param tolref=0.0\n');
fprintf(fid,'.param Tstop=%.6e\n.param tmax=0.25u\n',Tstop);
fclose(fid);

save('tx_info.mat','lI','lQ','tlead','Tsym','M','Fs','fc','N');
fprintf('Listo. Simula tx_lm1496.cir (Tstop = %.2f ms) y exporta V(sout) a salida.txt\n',Tstop*1e3);

function escribir_pwl(fn, v, tlead, Tsym)
    N  = numel(v); v = v(:);
    t0 = tlead + (0:N-1).'*Tsym;
    vprev = [v(1); v(1:end-1)];
    T = [0; reshape([t0.'; t0.'+20e-9],[],1); tlead+N*Tsym+1e-6];
    V = [v(1); reshape([vprev.'; v.'],[],1); v(end)];
    fid = fopen(fn,'w'); fprintf(fid,'%.9e %.6f\n',[T V].'); fclose(fid);
end
