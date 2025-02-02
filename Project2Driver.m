%% FCM Project 2
%Jonathan Engle
%Due Feb 20
%% Define Functions and parameters
clear;clc;clear all 
tic
%lower bound 
global a
a=-4;
%upper bound
global b
b=4;
%Number of steps
global n
n=5;
%Step size
global h
h=abs(a-b)/n;

%Normally spaced x vector
global nx
nx=linspace(a,b,n);
difx=linspace(a-1,b-1,n);

%Changing our active x makes it easier for testing
global activx
activx=nx; %chex%chexcheck; %chex lejaout %nx

%% Define functions
%Generic Function
fun=@(t)(t-2).^9;
funderiv=@(t) 9*(t+2).^8;%@(t) 9*t.^8-144*t.^7+ 1008*t.^6-4032*t.^5+10080*t.^4-16128*t.^3+16128*t.^2+9216*t;
%Runges Counter Example
run=@(t)(1./(1+t.^2));
runderiv=@(t)(-2*t/(1+t.^2).^2);
%fplot(run);

%% Cubic splines
[xi,yi]=parspline(activx,feval(run,activx));

error=norm(feval(run,activx))-norm(yi,inf);
%coefficients = cubic_spline(activx, feval(run,activx))
coefficients=yi;
hold on
plot(xi, yi)
fplot(run)
grid on
legend("Interpolated approximation","Test function")
hold off
%% Newtons basis
alphas = newton_interpolation(activx, feval(run,activx));
%Horners method
s=alphas;
    for i = n:1
        s(i) = s.*(xi-xi(i))+alphas(i);
    end

 error=norm(feval(fun,activx))-norm(s',inf);
feval(fun,activx);
% %Graphing
% hold on
% plot(activx,s);
% fplot(run);
% grid on
% legend("Interpolated approximation","Test function")
% hold off
%% hermite polynomial
[herm,hxv]=hermpol(activx,feval(fun,activx),feval(funderiv,activx),activx);
for i=1:n
    hermpoly(i)=-herm(i)*activx(i)^i;
end
hermpoly;
hxv;
%Graphing
% hold on
% plot(activx,herm);
% fplot(fun);
% grid on
% legend("Interpolated approximation","Test function")
% hold off
%  error=norm(feval(fun,difx)-herm,inf)%-norm(herm,inf)
toc