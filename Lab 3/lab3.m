%Ramaseshan Parthasarathy
%Lab 3

%%
%Problem 1 focuses on sinusoidal steady-state and transient response of
%filters.

%Problem 2 takes a sinsuoid input with some random noise, and processes it
%through a bandpass filter. Further more, zero-order hold and zero-mean
%white-noise signal sequence was investigated.

%%
%Problem 1
clear

u = @(t) (t>=0);
w0 = 4;

% part a
t_int = linspace(0,10,1001);
 
syms s t
x(t) = sin(w0*t)
H(s) = (s+3)/(s^2+s+1.25)
h(t) = ilaplace(H)
p = poles(H)
t_const = double(log(100)/abs(real(p(1))))

figure;
plot(t_int,h(t_int),'b');
grid on; xlabel('t(sec)');
title('impulse response, h(t)');
set(gca, 'XTick', 0:1:10, 'XLim', [0 10]);
set(gca, 'YTick', -2:1:2, 'YLim', [-2 2]);

%part b
eval = H(w0*1i)% H(jw)
mag_of_H = abs(H(w0 * 1i)).^2
phase = angle(eval)

yst = abs(H(w0*1i))*sin(w0*t_int + phase);

%part c
num = [1 3]; den = [1 1 1.25];
term = [1 -1i*w0]; den2 = conv(den, term)
[a,b,c] = residue(num, den2)

res = a(2) * exp(b(2)*t_int) + a(3) * exp(b(3)*t_int);
y_c = eval * exp(w0*1i*t_int) + res;
y = imag(y_c);

figure;
plot(t_int, x(t_int), 'k--', t_int, y, 'b');
grid on;
hold on; plot(t_int, yst, 'r--');
set(gca, 'XTick', 0:1:10, 'XLim', [0, 10]);
set(gca, 'YTick', -1:0.5:1, 'YLim', [-1, 1]);
title('y(t), y_{st}(t), x(t)');
xlabel('t(sec)');

%part d
tph = -phase/w0
x_val = x(t_int);
[t0, t_ind1] = max(x_val(800:900));
[t1, t_ind2] = max(y(800:900));

tmax_x = 800 + t_ind1 %x
tmax_y = 800 + t_ind2 %y
t_est = (tmax_y)/1000 - (tmax_x)/1000

%hold on; plot(tmax_x, y(tmax_x), 'r.', 'MarkerSize', 20);
%hold on; plot(tmax_y, y(tmax_y), 'r.', 'MarkerSize', 20);

%part e
ytr = y-yst;
figure; plot (t_int, ytr);
title ('transient, y_{tr}(t)');

clear all;
%%
%Problem 2

%part a
w0 = 5; alpha = 0.2;
w = 0:0.05:10;

H_mag_sq = @(w)(alpha^2 * w.^2)./((w.^2-w0^2).^2 + alpha^2 .* w.^2);

figure;
plot(w, H_mag_sq(w), 'b');
title('|H(j\omega)|^{2}, \omega_{0} = 5, \alpha = 0.2');
set(gca, 'XTick', 0:1:10, 'XLim', [0 10]);
set(gca, 'YTick', 0:0.5:1, 'YLim', [0 1.1]);
xlabel('\omega');

%part b
Tmax = 40; T = Tmax/2000;
t = 0:T:Tmax;
seed = 2016; rng(seed);
v = randn(size(t));
x = sin(w0*t) + v;

num = [alpha 0]; den = [1 alpha w0^2];
y = lsim(tf(num,den), x, t, [0;0], 'zoh');

figure; plot(t, x, 'k');
grid on; title('noisy input sinusoid, x(t)');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

figure; plot(t, y);
grid on; title('filtered output, y(t), \alpha = 0.2');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

%part c
y_v = lsim(tf(num,den),v,t,[0;0],'zoh');
figure; plot(t, v, 'k');
grid on; title('input noise, v(t)');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

figure; plot(t, y_v);
grid on; title('filtered noise, y_{v}(t), \alpha = 0.2');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

%part d

clear;

%part a
w0 = 5; alpha = 0.5;
w = 0:0.05:10

H_mag_sq = @(w)(alpha^2 * w.^2)./((w.^2-w0^2).^2 + alpha^2 .* w.^2);

figure;
plot(w, H_mag_sq(w), 'b');
title('|H(j\omega)|^{2}, \omega_{0} = 5, \alpha = 0.5');
set(gca, 'XTick', 0:1:10, 'XLim', [0 10]);
set(gca, 'YTick', 0:0.5:1, 'YLim', [0 1.1]);
xlabel('\omega');

%part b
Tmax = 40; T = Tmax/2000;
t = 0:T:Tmax;
seed = 2016; rng(seed);
v = randn(size(t));
x = sin(w0*t) + v;

num = [alpha 0]; den = [1 alpha w0^2];
y = lsim(tf(num,den), x, t, [0;0], 'zoh');

figure; plot(t, x, 'k');
grid on; title('noisy input sinusoid, x(t)');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

figure; plot(t, y);
grid on; title('filtered output, y(t), \alpha = 0.5');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

%part c
y_v = lsim(tf(num,den),v,t,[0;0],'zoh');
figure; plot(t, v, 'k');
grid on; title('input noise, v(t)');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

figure; plot(t, y_v);
grid on; title('filtered noise, y_{v}(t), \alpha = 0.5');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

clear;

%alpha = 1
%part a
w0 = 5; alpha = 1;
w = 0:0.05:10

H_mag_sq = @(w)(alpha^2 * w.^2)./((w.^2-w0^2).^2 + alpha^2 .* w.^2);

figure;
plot(w, H_mag_sq(w), 'b');
title('|H(j\omega)|^{2}, \omega_{0} = 5, \alpha = 1');
set(gca, 'XTick', 0:1:10, 'XLim', [0 10]);
set(gca, 'YTick', 0:0.5:1, 'YLim', [0 1.1]);
xlabel('\omega');

%part b
Tmax = 40; T = Tmax/2000;
t = 0:T:Tmax;
seed = 2016; rng(seed);
v = randn(size(t));
x = sin(w0*t) + v;

num = [alpha 0]; den = [1 alpha w0^2];
y = lsim(tf(num,den), x, t, [0;0], 'zoh');

figure; plot(t, x, 'k');
grid on; title('noisy input sinusoid, x(t)');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

figure; plot(t, y);
grid on; title('filtered output, y(t), \alpha = 1');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

%part c
y_v = lsim(tf(num,den),v,t,[0;0],'zoh');
figure; plot(t, v, 'k');
grid on; title('input noise, v(t)');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

figure; plot(t, y_v);
grid on; title('filtered noise, y_{v}(t), \alpha = 1');
set(gca, 'XTick', 0:10:40, 'XLim', [0, 40]);
set(gca, 'YTick', -4:1:4, 'YLim', [-4, 4]);

%part e
w0 = 5; alpha = 0.2;
Tmax = 40; T = Tmax/2000; t = 0:T:Tmax;
wr = sqrt((w0^2) - ((alpha^2)/4));

G = (alpha/wr) * exp(-alpha*(T/2)) * sin(wr*T);
a1 = -2 * exp(-alpha*(T/2)) * cos(wr*T); a2 = exp(-alpha*T);
Hd = @(z)(G.*z.^(-1) .* (1-z.^(-1)))./(1 + (a1.*z.^(-1)) + (a2.*z.^(-2)));
Hd(t);
xn = x;

ctr = 0; v1 = 0; v2 = 0; yn = 0;
while (ctr < t)
    yn = v1
    v1 = v2 + (G*xn) - (a1*yn);
    v2 = (-G*xn) - (a2*yn);
    ctr = ctr + 1;
end
disp(['yn = ' num2str(yn)]);

num = [0 G -G]; den = [1 a1 a2];
y_lsim = filter(num, den, x)
y_iter = yn;
E_lsim = norm(y_lsim - y_iter)
E_iter = norm(y_iter - y_lsim)

%part f
w0=5; alpha = 0.2; T = 0.02;
t_60dB = (log(1000)*2)/alpha;
N = (2 * t_60dB)/T;

gn0 = @(n)(alpha./wr) * exp((-alpha * n * T)./2) * sin(wr * n * T);
gn1 = @(n)(alpha./wr) * exp((-alpha * (n-1) * T)./2) * sin(wr * (n-1)* T);
hd = @(n)gn0(n) - gn1(n); p = exp((-alpha * T)/2) * exp(1i * wr *T);
op1 = G * exp((alpha * T)/2); op2 = (1 - exp((alpha*T)/2)) * exp(-1i * wr *T);
op3 = 2 * 1i * sin(wr * T); A = (op1 * op2)/op3;

NRR_exact = (2 * real((A.^2 * p.^2)./(1-(p.^2)))) + ...
            ((2 .* (abs(A).^2) * abs(p).^2)./(1-(abs(p).^2)))

ctr = 0; ans=0;
while (ctr < N-1)
    ans = (abs(hd(ctr)).^2) + ans;
    ctr = ctr + 1;
end

disp(['NRR1 = ' num2str(ans)])
NRR2 = (std(yn).^2)/(std(v).^2)
NRR3 = (T*alpha)/2
