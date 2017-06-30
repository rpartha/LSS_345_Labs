%Ramaseshan Parthasarathy
%Lab 2

%%
%Problem 1 serves to exhibit the following properties of LTI systems:
%linearity and time invariance. Further more, convolution was approximated
%numerically via the help of the conv method as used below.

%Problem 2 serves to exhbit transient and steady state (ss) responses of
%linear systems.

%%
%Problem 1

%part a
a = 0.9; T = 0.05; Tmax = 25;
t = [0:T:Tmax];  td = 1;
c1 = 1; c2 = 2; c3 = 1.5;
t1 = 0; t2 = 10; t3 = 15;
u = @(t) (t>=0);

h =  @(t) a.*exp(-a.*t).*u(t);
x = @(t,td) u(t) - u(t-td);


F = @(t, td) u(t) - u(t-td);
G = @(t,td) exp(-a*t) .* (exp(a*min(t,td)) - 1) .* (t>=0);

x1 =  c1*F(t-t1, td) + c2*F(t-t2, td) + c3*F(t-t3, td);
y1 =  c1*G(t-t1, td) + c2*G(t-t2, td) + c3*G(t-t3, td);

y  = T * conv(h(t), x1);
y = y(1:length(t));

figure;
plot(t, x1, 'k:', t, y, 'b-', t, y1, 'r--');
set(gca, 'Xtick', (0:5:25), 'XLim', [0 25]);
set(gca, 'YTick',(0:0.5:3), 'YLim', [0 3]);
title('t_{d} = 1, T = 0.005');
xlabel('t');
legend('input','conv', 'exact', 'Location','northeast');

clear;

%part b
a=0.9; T = 0.05; Tmax = 25;
t = [0:T:Tmax];  td = 3;
c1 = 1; c2 = 2; c3 = 1.5;
t1 = 0; t2 = 10; t3 = 15;
u = @(t) (t>=0);

h =  @(t) a.*exp(-a.*t).*u(t);
x = @(t,td) u(t) - u(t-td);


F = @(t, td) u(t) - u(t-td);
G = @(t,td) exp(-a*t) .* (exp(a*min(t,td)) - 1) .* (t>=0);

x1 =  c1*F(t-t1, td) + c2*F(t-t2, td) + c3*F(t-t3, td);
y1 =  c1*G(t-t1, td) + c2*G(t-t2, td) + c3*G(t-t3, td);

y  = T * conv(h(t), x1);
y = y(1:length(t));

figure;
plot(t, x1, 'k:', t, y, 'b-', t, y1, 'r--');
set(gca, 'Xtick', (0:5:25), 'XLim', [0 25]);
set(gca, 'YTick',(0:0.5:3), 'YLim', [0 3]);
title('t_{d} = 3, T = 0.005');
xlabel('t');
legend('input','conv', 'exact', 'Location','northeast');

clear;

a=0.9; T = 0.05; Tmax = 25;
t = [0:T:Tmax];  td = 5;
c1 = 1; c2 = 2; c3 = 1.5;
t1 = 0; t2 = 10; t3 = 15;
u = @(t) (t>=0);

h =  @(t) a.*exp(-a.*t).*u(t);
x = @(t,td) u(t) - u(t-td);


F = @(t, td) u(t) - u(t-td);
G = @(t,td) exp(-a*t) .* (exp(a*min(t,td)) - 1) .* (t>=0);

x1 =  c1*F(t-t1, td) + c2*F(t-t2, td) + c3*F(t-t3, td);
y1 =  c1*G(t-t1, td) + c2*G(t-t2, td) + c3*G(t-t3, td);

y  = T * conv(h(t), x1);
y = y(1:length(t));

figure;
plot(t, x1, 'k:', t, y, 'b-', t, y1, 'r--');
set(gca, 'Xtick', (0:5:25), 'XLim', [0 25]);
set(gca, 'YTick',(0:0.5:3), 'YLim', [0 3]);
title('t_{d} = 5, T = 0.005');
xlabel('t');
legend('input','conv', 'exact', 'Location','northeast');

clear;

%part c

%part a repeated
a=0.9; T = 0.01; Tmax = 25;
t = [0:T:Tmax];  td = 1;
c1 = 1; c2 = 2; c3 = 1.5;
t1 = 0; t2 = 10; t3 = 15;
u = @(t) (t>=0);

h =  @(t) a.*exp(-a.*t).*u(t);
x = @(t,td) u(t) - u(t-td);


F = @(t, td) u(t) - u(t-td);
G = @(t,td) exp(-a*t) .* (exp(a*min(t,td)) - 1) .* (t>=0);

x1 =  c1*F(t-t1, td) + c2*F(t-t2, td) + c3*F(t-t3, td);
y1 =  c1*G(t-t1, td) + c2*G(t-t2, td) + c3*G(t-t3, td);

y  = T * conv(h(t), x1);
y = y(1:length(t));

figure;
plot(t, x1, 'k:', t, y, 'b-', t, y1, 'r--');
set(gca, 'Xtick', (0:5:25), 'XLim', [0 25]);
set(gca, 'YTick',(0:0.5:3), 'YLim', [0 3]);
title('t_{d} = 1, T = 0.001');
xlabel('t');
legend('input','conv', 'exact', 'Location','northeast');

clear;

%part b repeated
a=0.9; T = 0.01; Tmax = 25;
t = [0:T:Tmax];  td = 3;
c1 = 1; c2 = 2; c3 = 1.5;
t1 = 0; t2 = 10; t3 = 15;
u = @(t) (t>=0);

h =  @(t) a.*exp(-a.*t).*u(t);
x = @(t,td) u(t) - u(t-td);


F = @(t, td) u(t) - u(t-td);
G = @(t,td) exp(-a*t) .* (exp(a*min(t,td)) - 1) .* (t>=0);

x1 =  c1*F(t-t1, td) + c2*F(t-t2, td) + c3*F(t-t3, td);
y1 =  c1*G(t-t1, td) + c2*G(t-t2, td) + c3*G(t-t3, td);

y  = T * conv(h(t), x1);
y = y(1:length(t));

figure;
plot(t, x1, 'k:', t, y, 'b-', t, y1, 'r--');
set(gca, 'Xtick', (0:5:25), 'XLim', [0 25]);
set(gca, 'YTick',(0:0.5:3), 'YLim', [0 3]);
title('t_{d} = 3, T = 0.001');
xlabel('t');
legend('input','conv', 'exact', 'Location','northeast');

clear;

a=0.9; T = 0.01; Tmax = 25;
t = [0:T:Tmax];  td = 5;
c1 = 1; c2 = 2; c3 = 1.5;
t1 = 0; t2 = 10; t3 = 15;
u = @(t) (t>=0);

h =  @(t) a.*exp(-a.*t).*u(t);
x = @(t,td) u(t) - u(t-td);


F = @(t, td) u(t) - u(t-td);
G = @(t,td) exp(-a*t) .* (exp(a*min(t,td)) - 1) .* (t>=0);

x1 =  c1*F(t-t1, td) + c2*F(t-t2, td) + c3*F(t-t3, td);
y1 =  c1*G(t-t1, td) + c2*G(t-t2, td) + c3*G(t-t3, td);

y  = T * conv(h(t), x1);
y = y(1:length(t));

figure;
plot(t, x1, 'k:', t, y, 'b-', t, y1, 'r--');
set(gca, 'Xtick', (0:5:25), 'XLim', [0 25]);
set(gca, 'YTick',(0:0.5:3), 'YLim', [0 3]);
title('t_{d} = 5, T = 0.001');
xlabel('t');
legend('input','conv', 'exact', 'Location','northeast');

clear all;

%%
%Problem 2

%part a
w0 = 2; w1 = 3; alpha = 0.3;
u = @(t) (t>=0);
F = @(t, td) u(t) - u(t-td);
Tmax = 100; T = Tmax/2000;
t = 0:T:Tmax;


x = sin(w1*t) .* F(t, 30) + ...
    sin(w0*t) .* F(t-30, 40) + ...
    sin(w1*t) .* F(t-70, 30);

wr = sqrt((w0^2) - ((alpha^2)/4));
g = (alpha*exp((-alpha*t)/2)).*(cos(wr*t) - (alpha/(2*wr))*sin(wr*t)).*u(t);
con = conv(g,x);
y = x - T * con(1:length(t));

H = @(s) (s.^2 + w0^2)./(s.^2 + alpha*s + w0^2);
val = abs(H(1i * w1));

figure;
plot(t, x, 'b');
set(gca, 'Ytick', (-2:1:2), 'YLim', [-2 2]);
set(gca, 'XTick', (0:10:100), 'XLim', [0 100]);
title('input signal, x(t)'); grid on
xlabel('t');

figure;
t1 = 0:T:30; t2 = 70:T:100;
y1 = val * ones(1, length(t1));
y2 = val * ones(1, length(t2));

plot(t, y, 'r');
set(gca, 'Ytick', (-2:1:2), 'YLim', [-2 2]);
set(gca, 'XTick', (0:10:100), 'XLim', [0 100]);
title('output signal, y(t), conv method, \alpha = 0.3'); grid on
xlabel('t');
hold on;
plot(t1, y1,'b', t2, y2, 'b', t1, -y1, 'b', t2, -y2, 'b');
legend('y(t)', '|H(jw_{1})|', 'Location', 'north');

%part b
p = (-alpha/2) + 1i*wr;
tau = log(100)/abs(real(p))

%part c
num = [1, 0, w0^2]; den = [1, alpha, w0^2];
tfn = tf(num, den);
y4 = (lsim(tfn, x, t))';

figure;
plot(t, y4,'r-', t1, y1,'b-', t2, y2, 'b-', t1, -y1, 'b-', t2, -y2, 'b-');
set(gca, 'Ytick', (-2:1:2), 'YLim', [-2 2]);
set(gca, 'XTick', (0:10:100), 'XLim', [0 100]);
title('output signal, y(t), lsim method, \alpha = 0.3'); grid on
legend('y(t)', '|H(jw_{1})|', 'Location', 'north');
xlabel('t');

p_err = 100 * (norm(y - y4)/norm(y)) %percent error

%part d
one = w0^2 + (alpha^2)/2;
two = (alpha * sqrt(w0^2)) + (alpha^2)/4;

w_plus = sqrt(one + two);
w_minus = sqrt(one - two);

alpha_def = w_plus - w_minus;
w0_squared = w_plus * w_minus;

w_int = 0:0.01:5;
mag_of_H = abs(H(w_int * 1i)).^2; % H(jw)

three_dB_frequency = w_minus:0.01:w_plus;
mag_half_order = 0.5 * ones(1, length(three_dB_frequency));

figure;
plot(w_int, mag_of_H, three_dB_frequency, ...
     mag_half_order, w1, mag_of_H(find(w_int == w1)), ...
     'r.', 'MarkerSize', 20);
hold on; plot(2, 0, 'r.', 'MarkerSize', 20);
set(gca, 'Ytick', (0:0.5:1), 'YLim', [0 1.15]);
set(gca, 'XTick', (0:1:5), 'XLim', [0 5]);
title('notch filter, |H(\omega)|^{2}, \alpha = 0.3');
text(2.5, 0.5, '\color{red} 3-dB width');
text(w1-0.085, mag_of_H(find(w_int == w1))-0.05, '\bf \uparrow');
text(w1-0.085, mag_of_H(find(w_int == w1))-0.1, '\bf \omega_{1}');

clear all;
%parts a-d repeated
%part a
w0 = 2; w1 = 3; alpha = 0.9;
u = @(t) (t>=0);
F = @(t, td) u(t) - u(t-td);
Tmax = 100; T = Tmax/2000;
t = 0:T:Tmax;


x = sin(w1*t) .* F(t, 30) + ...
    sin(w0*t) .* F(t-30, 40) + ...
    sin(w1*t) .* F(t-70, 30);

wr = sqrt((w0^2) - ((alpha^2)/4));
g = (alpha*exp((-alpha*t)/2)).*(cos(wr*t) - (alpha/(2*wr))*sin(wr*t)).*u(t);
con = conv(g,x);
y = x - T * con(1:length(t));

H = @(s) (s.^2 + w0^2)./(s.^2 + alpha*s + w0^2);
val = abs(H(1i * w1));

figure;
plot(t, x, 'b');
set(gca, 'Ytick', (-2:1:2), 'YLim', [-2 2]);
set(gca, 'XTick', (0:10:100), 'XLim', [0 100]);
title('input signal, x(t)'); grid on
xlabel('t');

figure;
t1 = 0:T:30; t2 = 70:T:100;
y1 = val * ones(1, length(t1));
y2 = val * ones(1, length(t2));

plot(t, y, 'r');
set(gca, 'Ytick', (-2:1:2), 'YLim', [-2 2]);
set(gca, 'XTick', (0:10:100), 'XLim', [0 100]);
title('output signal, y(t), conv method, \alpha = 0.9'); grid on
xlabel('t');
hold on;
plot(t1, y1,'b', t2, y2, 'b', t1, -y1, 'b', t2, -y2, 'b');
legend('y(t)', '|H(jw_{1})|', 'Location', 'north');

%part b
p = (-alpha/2) + 1i*wr;
tau = log(100)/abs(real(p))

%part c
num = [1, 0, w0^2]; den = [1, alpha, w0^2];
tfn = tf(num, den);
y4 = (lsim(tfn, x, t))';

figure;
plot(t, y4,'r-', t1, y1,'b-', t2, y2, 'b-', t1, -y1, 'b-', t2, -y2, 'b-');
set(gca, 'Ytick', (-2:1:2), 'YLim', [-2 2]);
set(gca, 'XTick', (0:10:100), 'XLim', [0 100]);
title('output signal, y(t), lsim method, \alpha = 0.9'); grid on
legend('y(t)', '|H(jw_{1})|', 'Location', 'north');
xlabel('t');

p_err = 100 * (norm(y - y4)/norm(y)) %percent error

%part d
one = w0^2 + (alpha^2)/2;
two = (alpha * sqrt(w0^2)) + (alpha^2)/4;

w_plus = sqrt(one + two);
w_minus = sqrt(one - two);

alpha_def = w_plus - w_minus;
w0_squared = w_plus * w_minus;

w_int = 0:0.01:5;
mag_of_H = abs(H(w_int * 1i)).^2; % H(jw)

three_dB_frequency = w_minus:0.01:w_plus;
mag_half_order = 0.5 * ones(1, length(three_dB_frequency));

figure;
plot(w_int, mag_of_H, three_dB_frequency, ...
     mag_half_order, w1, mag_of_H(find(w_int == w1)), ...
     'r.', 'MarkerSize', 20);
hold on; plot(2, 0, 'r.', 'MarkerSize', 20);
set(gca, 'Ytick', (0:0.5:1), 'YLim', [0 1.15]);
set(gca, 'XTick', (0:1:5), 'XLim', [0 5]);
title('notch filter, |H(\omega)|^{2}, \alpha = 0.9');
text(2.7, 0.5, '\color{red} 3-dB width');
text(w1-0.085, mag_of_H(find(w_int == w1))-0.05, '\bf \uparrow');
text(w1-0.085, mag_of_H(find(w_int == w1))-0.1, '\bf \omega_{1}');

clear all;
%%
%Problem 3

Y = load('lab2.dat');

%part a
f0 = 60; width_3dB = 1.5;
w0 = 2 * pi * f0;
alpha = 2 * pi * width_3dB;

one = w0^2 + (alpha^2)/2;
two = (alpha * sqrt(w0^2)) + (alpha^2)/4;
w_left = sqrt(one + two)
w_right = sqrt(one - two)

Q = f0/width_3dB %quality factor

wr = sqrt((w0^2) - ((alpha^2)/4));
p = (-alpha/2) + 1i*wr;
tau = log(100)/abs(real(p))

%part b
t = Y(:,1);
s_t = Y(:,2) ;
x_t = Y(:,3);

figure; plot(t, s_t, 'k');
set(gca, 'Ytick', (-1:0.5:1), 'YLim', [-1 1.15]);
set(gca, 'XTick', (0:0.5:1.5), 'XLim', [0 1.5]);
title('noise-free ECG, s(t)');
xlabel('t (sec)');

figure; plot(t, x_t);
set(gca, 'Ytick', (-1:0.5:1), 'YLim', [-1.15 1.15]);
set(gca, 'XTick', (0:0.5:1.5), 'XLim', [0 1.5]);
title('ECG + 60-Hz interference, x(t)');
xlabel('t (sec)');

%part c
num = [1, 0, w0^2]; den = [1, alpha, w0^2];
tfn = tf(num,den);

y = (lsim(tfn, x_t, t))';
figure; plot(t, y, 'r', t, s_t, 'k');
set(gca, 'Ytick', (-1:0.5:1), 'YLim', [-1.15 1.15]);
set(gca, 'XTick', (0:0.5:1.5), 'XLim', [0 1.5]);
legend('y(t)', 's(t)', 'Location', 'southeast');
title('processed ECG, y(t)');
xlabel('t (sec)');

%part d
w_int = 0:1:120;
H = @(s) (s.^2 + w0^2)./(s.^2 + alpha*s + w0^2);
mag_of_H = abs(H(w_int * 2 * pi * 1i)).^2;

figure;
plot(w_int, mag_of_H, 'b');
set(gca, 'Ytick', (0:0.5:1), 'YLim', [0 1.15]);
set(gca, 'XTick', (0:30:120), 'XLim', [0 120]);
title('notch-filter magnitude response, |H(f)|^{2}');
xlabel('t (sec)');
