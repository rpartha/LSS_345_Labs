%Ramaseshan Parthasarathy
%September 26, 2016

%LAB 1
%PURPOSE AND METHODS USED:
%For part (1) of this lab, four limiting forms of the dirac-delta function
%were plotted to depict its definition. This was accomplished
%by implementing anonymous functions.
%For the second part of this lab, a linear system was given with an intital
%condition y(0-), and was replaced by a state-space model which was the
%same. The following were accomplished:
 %1. The impulse response is found by using the method of
 %ilaplace function, which yields h(t).
 %2. Knowing the non-zero initial condition, y(t) is found by lapacing,
 %partial fraction expansion, and inversing it. Repeat for v(t).
 %3. Alternatively solve for y(t) and v(t) the dsolve method.
 %4. However, apart from analytical methods previously described, y(t)
 %and v(t) were both found numerically using the lsim method.
 %This requires the tf to passed in a state-space objects.
 %5. Lastly, we consider the discrete-time implementation of the system.
 %6. Processes 1-5 are again repeated for a initial condition of zero.

%PROBLEM 1
t = linspace(-1, 1, 1001);
u = @(t) (t>=0);
d1 = @(t,e) (1/e)*(u(t + 0.5*e) - u(t - 0.5*e));
figure; plot(t, d1(t,.05), t, d1(t, .10), 'r:');
set(gca, 'YTick', 0:5:25, 'YLim', [0 25]);
set(gca, 'XTick', -1:.25:1, 'XLim', [-1 1]);
title('rectangular pulse');
xlabel('t');
legend('\epsilon = 0.05','\epsilon = 0.10','Location','northeast');

d2= @(t,e) (1/sqrt(2 * pi * e)) * exp(-(t.^2)/(2*e));
figure; plot(t, d2(t,.001), t, d2(t, .002), 'r:');
set(gca, 'YTick', 0:5:15, 'YLim', [0 15]);
set(gca, 'XTick', -1:.25:1, 'XLim', [-1 1]);
title('gaussian');
xlabel('t');
legend('\epsilon = 0.001','\epsilon = 0.002','Location','northeast');

d3 = @(t,e) (1/pi) * (e./((e^2) + (t.^2)));
figure; plot(t, d3(t,.005), t, d3(t, .010), 'r:');
set(gca, 'YTick', 0:5:70, 'YLim', [0 70]);
set(gca, 'XTick', -1:.25:1, 'XLim', [-1 1]);
title('\delta_{3}(t)');
xlabel('t');
legend('\epsilon = 0.005','\epsilon = 0.010','Location','northeast');

d4= @(t,e) sinc(t/e)./(pi*t);
figure; plot(t, d4(t,.05), t, d4(t, .010), 'r:');
set(gca, 'YTick', -40:5:50, 'YLim', [-40 50]);
set(gca, 'XTick', -1:.25:1, 'XLim', [-1 1]);
title('\delta_{4}(t)');
xlabel('t');
legend('\epsilon = 0.05','\epsilon = 0.010','Location','northeast');
clear all;

%PROBLEM 2
%close all;
clear all;

%part a
%--see part (f)--%
syms s t Y V;

%part b
H = (3*s + 2)/(s + 2);
h = ilaplace(H) %3*dirac(t)- 4e^(-2t)

%part c
v0 = -1/-4; %v(0+)= v(0-) = v0 = y(0-)/(b1 - a*b0)= -1/(2 - 3*2)

%part d
y0 = -1; b1 = 2; a = 2; b0 = 3;
x = exp(-t);
x0 = subs(x,t,0);
X = laplace(x);
Y = solve(s*Y-y0 + a*Y == b0*s*X + b1*X, Y); %Y = (2s + 1)/(s^2 + 3s + 2)
Y = partfrac(Y, s); %Y = 3/(s+2) - 1/(s+1)
y1 = ilaplace(Y) %y = 3e^(-2t) - e^(-t)

%part e
syms t y(t)
dy = diff(y,t); dx = diff(x,t);
y2 = dsolve(dy + a*y == b0*dx + b1*x, y(0) == y0 + b0*x0) %method 2
simplify(y1-y2); %check

%part f
V = solve(s*V - v0 + 2*V == X, V);
V = partfrac(V);
v1 = ilaplace(V)
syms t v(t)
dv = diff(v,t); dx = diff(x,t);
%v2 = dsolve(dv + 2*v1 == b0*dx + b1*x, v(0) == v0 + b0*x0)
v2 = dsolve(dv + 2*v == x, v(0) == v0); %method 2
simplify(y1 + (4*v1) - (3*x)); %check ans = 0

%part g
t1 = linspace(0, 5, 201);
y3 = subs(y1, symvar(y1), t1);
v3 = subs(v1,symvar(v1), t1);
figure; plot(t1, y3, 'b-', t1, v3, 'r--'); %plot exact expressions
set(gca, 'XTick', 0:1:5, 'XLim', [0 5]);
%set(gca, 'YTick', 0:.25:2, 'YLim', [0 2]);
title('exact outputs, y(0^{-}) = -1');
xlabel('t');
legend('y(t)','v(t)','Location','northeast');
clear;

%part h
num = [3,2]; den = [1,2];
H = tf(num,den);
[a,b,c,d] = tf2ss(num,den);
S = ss(a,b,c,d);
t2 = linspace(0,5,201)';
x = exp(-t2);
s = tf('s');
Hv = 1/(s+1);
Sv = ss(Hv);
v0_1 = -1/-4;
y4 = lsim(S,x,t2,v0_1);
v4 = lsim(Sv,x,t2,v0_1);
norm(y4 + 4*v4 - 3*x)
figure; plot(t2, y4, 'b-', t2, v4, 'r--'); %plot exact expressions
set(gca, 'XTick', 0:1:5, 'XLim', [0 5]);
%set(gca, 'YTick', 0:.25:2, 'YLim', [0 2]);
title('lsim outputs, y(0^{-}) = -1');
xlabel('t');
legend('y(t)','v(t)','Location','northeast');
clear all;

%part i
%discrete-time method
T = 0.1;
tn = 0:T:5;
N = length(tn);
a = 2; b0 = 3; b1 = 2;
y0 = -1; v0 = y0/(b1-a*b0);
a1 = -exp(-a*T);
b0_new = b0;
b1_new = b1*(1-exp(-a*T))/a - b0_new
x = exp(-2*tn);
w = y0; v = 0;
for n=0:N-1,
 y(n+1) = -a1*w + b0_new*x(n+1) + b1_new*v;
 w = y(n+1);
 v = x(n+1);
end
w = 0; v = 0;
clear v
a1 = -exp(-a*T);
b1_new = (1 - exp(-a*T))/a;
w = v0; u = 0;
for n=0:N-1,
 v(n+1) = -a1*w + b1_new*u;
 w = v(n+1);
 u = x(n+1);
end
%exact values obtained from part (g)
syms t
y1 = 3*exp(-2*t)-exp(-t);
v1 = exp(-t) - (3*exp(-2*t))/4;
t1 = linspace(0,5,201);
y1a = subs(y1, t, t1);
v1a = subs(v1, t, t1);
figure;
plot(tn, y, 'b-', tn, v, 'g-', t1, y1a, 'r:', t1,
 v1a, 'r:', 'linewidth', 1);
ylim([-.25, 2.25]);
7
xlim([0,5]);
set(gca, 'XTick', [0 1 2 3 4 5]);
set(gca, 'YTick', [0 0.5 1 1.5 2]);
title('discrete-time outputs, y(0^{-}) = -1, T = 0.1');
legend('v(t_{n})', 'y(t_{n})', 'exact', 'Location', 'northeast');
clear all;

%part j
%--repeat of parts (c)-(i) when y(0-) = 0--%
%part c
syms s t Y V;
v0 = 0; %v(0+)= v(0-) = v0 = y(0-)/(b1 - a*b0)= 0/(2 - 3*2)
%part d
y0 = 0; b1 = 2; a = 2; b0 = 3;
x = exp(-t);
x0 = subs(x,t,0);
X = laplace(x);
Y = solve(s*Y-y0 + a*Y == b0*s*X + b1*X, Y); %Y = (2s + 1)/(s^2 + 3s + 2)
Y = partfrac(Y, s); %Y = 3/(s+2) - 1/(s+1)
y1 = ilaplace(Y) %y = 3e^(-2t) - e^(-t)
%part e
syms t y(t)
dy = diff(y,t); dx = diff(x,t);
y2 = dsolve(dy + a*y == b0*dx + b1*x, y(0) == y0 + b0*x0) %method 2
simplify(y1-y2); %check
%part f
V = solve(s*V - v0 + 2*V == X, V);
V = partfrac(V);
v1 = ilaplace(V)
syms t v(t)
dv = diff(v,t); dx = diff(x,t);
%v2 = dsolve(dv + 2*v1 == b0*dx + b1*x, v(0) == v0 + b0*x0)
v2 = dsolve(dv + 2*v == x, v(0) == v0); %method 2
simplify(y1 + (4*v1) - (3*x)); %check ans = 0
%part g
t1 = linspace(0, 5, 201);
y3 = subs(y1, symvar(y1), t1);
v3 = subs(v1,symvar(v1), t1);
figure; plot(t1, y3, 'b-', t1, v3, 'r--'); %plot exact expressions
set(gca, 'XTick', 0:1:5, 'XLim', [0 5]);
%set(gca, 'YTick', 0:.25:2, 'YLim', [0 2]);
title('exact outputs, y(0^{-}) = -1');
xlabel('t');
legend('y(t)','v(t)','Location','northeast');
clear;
%part h
num = [3,2]; den = [1,2];
H = tf(num,den);
[a,b,c,d] = tf2ss(num,den);
S = ss(a,b,c,d);
t2 = linspace(0,5,201)';
x = exp(-t2);
s = tf('s');
Hv = 1/(s+1);
Sv = ss(Hv);
v0_1 = 0;
y4 = lsim(S, x, t2, v0_1);
v4 = lsim(Sv, x, t2, v0_1);
norm(y4 + 4*v4 - 3*x)
figure; plot(t2, y4, 'b-', t2, v4, 'r--'); %plot exact expressions
set(gca, 'XTick', 0:1:5, 'XLim', [0 5]);
%set(gca, 'YTick', 0:.25:2, 'YLim', [0 2]);
title('lsim outputs, y(0^{-}) = -1');
xlabel('t');
legend('y(t)','v(t)','Location','northeast');
clear;
%part i
%discrete-time method
T = 0.1;
tn = 0:T:5;
N = length(tn);
a = 2; b0 = 3; b1 = 2;
y0 = 0; v0 = y0/(b1-a*b0);
a1 = -exp(-a*T);
b0_new = b0;
b1_new = b1*(1-exp(-a*T))/a - b0_new
x = exp(-2*tn);
w = y0; v = 0;
for n=0:N-1,
 y(n+1) = -a1*w + b0_new*x(n+1) + b1_new*v;
 w = y(n+1);
 v = x(n+1);
end
w = 0; v = 0;
clear v
a1 = -exp(-a*T);
b1_new = (1 - exp(-a*T))/a;
w = v0; u = 0;
for n=0:N-1,
 v(n+1) = -a1*w + b1_new*u;
 w = v(n+1);
 u = x(n+1);
end
%exact values obtained from part (g)
syms t
v1 = exp(-t) - exp(-2*t);
y1 = 4*exp(-2*t) - exp(-t);
t1 = linspace(0,5,201);
y1a = subs(y1, t, t1);
v1a = subs(v1, t, t1);
figure;
plot(tn, y, 'b-', tn, v, 'g-', t1, y1a, 'r:', t1,
 v1a, 'r:', 'linewidth', 1);
ylim([-.25, 2.25]);
xlim([0,5]);
set(gca, 'XTick', [0 1 2 3 4 5]);
set(gca, 'YTick', [0 0.5 1 1.5 2]);
title('discrete-time outputs, y(0^{-}) = 0, T = 0.1');
legend('v(t_{n})', 'y(t_{n})', 'exact', 'Location', 'northeast');
