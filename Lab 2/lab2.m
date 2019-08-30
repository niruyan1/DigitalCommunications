%% Niruyan Rakulan 214343438 EECS 4214 Lab 2
%% Problem 1
clc;
close all;
clear all;

%Time Domain
t=0:1e-6:1;
f1=1e3;
x=sin(2*pi*f1*t)+(1/3)*sin(2*pi*3*f1*t)+(1/5)*sin(2*pi*5*f1*t);
subplot(2,1,1);
plot(t,x);
axis([0 0.999e-3 -1.5 1.5]);
title('x(t)');
xlabel('Time (s)');
ylabel('Amp.');
grid; 
%Freq. Domain
subplot(2,1,2);

%Normalize FFT
X=fft(x);
X=X/length(x);
X = X(1:length(X)/2+1);
X(2:end-1) = 2*X(2:end-1);
f = linspace(0,1e6/2,length(X));
plot(f,abs(X));
axis([0 10000 0 2]);
title('X(jw)');
xlabel('Freq.(Hz)');
ylabel('Amp.');
grid;

%% Problem 2
clc;
close all;
clear all;

%Time Domain
t=0:1e-6:1;
f1=1e3;
x=sin(2*pi*f1*t)+(1/3)*sin(2*pi*3*f1*t)+(1/5)*sin(2*pi*5*f1*t);
x=x+2*randn(1,length(t));
subplot(2,1,1);
plot(t,x);
axis([0 0.999e-3 -10 10]);
title('x(t)');
xlabel('Time (s)');
ylabel('Amp.');

%Freq. Domain
subplot(2,1,2);

X=fft(x);
X=X/length(x);
X = X(1:length(X)/2+1);
X(2:end-1) = 2*X(2:end-1);
f = linspace(0,1e6/2,length(X));
plot(f,abs(X));
axis([0 10000 0 2]);
title('X(jw)');
xlabel('Freq.(Hz)');
ylabel('Amp.');
%% Problem 3
clc;
close all;
clear all;
vars=randn(100,1000);

%time averages
ta=mean(vars(ceil(rand*100),:));
tb=mean(vars(ceil(rand*100),:));
tc=mean(vars(ceil(rand*100),:));

%ensemble averages
ea=mean(vars(:,ceil(rand*1000)));
eb=mean(vars(:,ceil(rand*1000)));
ec=mean(vars(:,ceil(rand*1000)));

fprintf('Time Average 1: %f\n', ta);
fprintf('Time Average 2: %f\n', tb);
fprintf('Time Average 3: %f\n', tc);

fprintf('Ensemble Average 1: %f\n', ea);
fprintf('Ensemble Average 1: %f\n', eb);
fprintf('Ensemble Average 1: %f\n', ec);

fprintf('The ensemble values are fairly close to the time averages (near 0) considering the sample size; therefore the process is ergodic. \n');

%% Problem 4
clc;
close all;
clear all;

%animation can skip
figure();
%%WARNING: CTRL_C to quit out of animation
syms y(t);
y(t) = piecewise(1<t<3, 1, 3<t<4, -1, 0);
subplot(2,1,1);
fplot(y(t),'r');
hold on;
axis([-4.5 9.5 -2 2]);
ylabel('Amplitude');
xlabel('Time(s)');
title('x(t)');
subplot(2,1,2);
hold on;
axis([-5 5 -2 4]);
ylabel('Amplitude');
xlabel('T(s)');
title('R(T)');
for shift=-5:0.1:5
    subplot(2,1,1);
    f=fplot(y(t+shift),'b');
    legend('x(t)','x(t+T)');
    subplot(2,1,2);
    sumz=0;
    for count=-4.5:0.1:9.5
        if (y(count)~=0&&y(count+shift)~=0)
            sumz=(y(count)*y(count+shift))*0.1+sumz;
        end
    end
    plot(shift,sumz,'g.');
    pause(0.01);
    delete(f);
end

%Rx plotted using piecewise function
figure();
syms Rx(T);
Rx(T) = piecewise(-inf<T<-3, 0, -3<T<-2, -3-T, -2<T<-1, 1+T, -1<T<0, 3+3*T,0<=T<1, 3-3*T, 1<T<2, 1-T, 2<T<3, -3+T, 0);
fplot(Rx(T));
axis([-5 5 -2 4]);
ylabel('Amplitude');
xlabel('T(s)');
title('R(T)');

%Hand written
figure();
imshow('autocorrelation.png');
