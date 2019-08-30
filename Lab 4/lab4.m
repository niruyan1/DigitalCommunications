%% Niruyan Rakulan 214343438 Lab 4
%% Problem 1
%1.1
clc;
close all;
clear all;

t=0:1e-9:(2e-3)-(1e-9);
x=5*sawtooth(2*pi*1e3*t);
quantile_interval=10/(2^4);
quantizer=floor(x/quantile_interval)*quantile_interval;
figure;
subplot(2,1,1);
plot(t,x);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Orignal Signal');
grid;
axis([0 2e-3 -6 6]);

subplot(2,1,2);
plot(t,quantizer);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Quantized Signal (4 bits)');
grid;
axis([0 2e-3 -6 6]);

%1.2
figure;
error=x-quantizer;
plot(t,error);
xlabel('Time(s)');
ylabel('Error Voltage(V)');
title('Error (4 bits)');
axis([0 2e-3 -0.75 0.75]);

%1.3 RMS is the same as last method
rms=sqrt(mean(error.^2));

fprintf('1.3: RMS Value= %f\n',rms);

%% QEO 4bits
clc;
close all;
clear all;

t=0:1e-9:(2e-3)-(1e-9);
x=5*sawtooth(2*pi*1e3*t);
quantile_interval=10/(2^4);
quantizer=floor(x/quantile_interval)*quantile_interval+(quantile_interval/2);
figure;
subplot(2,1,1);
plot(t,x);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Orignal Signal');
grid;
axis([0 2e-3 -6 6]);

subplot(2,1,2);
plot(t,quantizer);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Quantized Signal (4 bits QEO)');
grid;
axis([0 2e-3 -6 6]);

figure;
error=x-quantizer;
plot(t,error);
xlabel('Time(s)');
ylabel('Error Voltage(V)');
title('Error (4 bits QEO)');
axis([0 2e-3 -0.75 0.75]);
rms=sqrt(mean(error.^2));



quantization_noise_power=(quantile_interval)^2/12;

fprintf('4 bits QEO: RMS Value= %f\nQuantization Noise Power= %f.\nThe RMS value is half the value of the old method(0.360844)\n',rms, quantization_noise_power);

figure;
histogram(error);
xlabel('Error(V)');
ylabel('Occurences');
title('Histogram of Errors');
%% 6 bits
clc;
close all;
clear all;

t=0:1e-9:(2e-3)-(1e-9);
x=5*sawtooth(2*pi*1e3*t);
quantile_interval=10/(2^6);
quantizer=floor(x/quantile_interval)*quantile_interval+(quantile_interval/2);

figure;
subplot(2,1,1);
plot(t,x);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Orignal Signal');
grid;
axis([0 2e-3 -6 6]);

subplot(2,1,2);
plot(t,quantizer);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Quantized Signal(6 bits QEO)');
grid;
axis([0 2e-3 -6 6]);

figure;
error=x-quantizer;
plot(t,error);
xlabel('Time(s)');
ylabel('Error Voltage(V)');
title('Error (6 bits QEO)');
axis([0 2e-3 -0.75 0.75]);
rms=sqrt(mean(error.^2));

quantization_noise_power=(quantile_interval)^2/12;

fprintf('6 bits QEO: RMS Value= %f\nQuantization Noise Power= %f\n',rms, quantization_noise_power);

figure;
histogram(error);
xlabel('Error(V)');
ylabel('Occurences');
title('Histogram of Errors');
%% 8 bits
clc;
close all;
clear all;

t=0:1e-9:(2e-3)-(1e-9);
x=5*sawtooth(2*pi*1e3*t);
quantile_interval=10/(2^8);
quantizer=floor(x/quantile_interval)*quantile_interval+(quantile_interval/2);

figure;
subplot(2,1,1);
plot(t,x);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Orignal Signal');
grid;
axis([0 2e-3 -6 6]);

subplot(2,1,2);
plot(t,quantizer);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Quantized Signal(8 bits QEO)');
grid;
axis([0 2e-3 -6 6]);

figure;
error=x-quantizer;
plot(t,error);
xlabel('Time(s)');
ylabel('Error Voltage(V)');
title('Error (8 bits QEO)');
axis([0 2e-3 -0.75 0.75]);
rms=sqrt(mean(error.^2));

quantization_noise_power=(quantile_interval)^2/12;

fprintf('8 bits QEO: RMS Value= %f\nQuantization Noise Power= %f\n',rms, quantization_noise_power);

figure;
histogram(error);
xlabel('Error(V)');
ylabel('Occurences');
title('Histogram of Errors');
%% 10 bits
clc;
close all;
clear all;

t=0:1e-9:(2e-3)-(1e-9);
x=5*sawtooth(2*pi*1e3*t);
quantile_interval=10/(2^10);
quantizer=floor(x/quantile_interval)*quantile_interval+(quantile_interval/2);

figure;
subplot(2,1,1);
plot(t,x);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Orignal Signal');
grid;
axis([0 2e-3 -6 6]);

subplot(2,1,2);
plot(t,quantizer);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Quantized Signal(10 bits QEO)');
grid;
axis([0 2e-3 -6 6]);

figure;
error=x-quantizer;
plot(t,error);
xlabel('Time(s)');
ylabel('Error Voltage(V)');
title('Error(10 bits QEO)');
axis([0 2e-3 -0.75 0.75]);
rms=sqrt(mean(error.^2));

quantization_noise_power=(quantile_interval)^2/12;

fprintf('10 bits QEO: RMS Value= %f\nQuantization Noise Power= %f\n',rms, quantization_noise_power);

figure;
histogram(error);
xlabel('Error(V)');
ylabel('Occurences');
title('Histogram of Errors');

%% 1.6
clc;
close all;
clear all;

disp('4 bits QEO: RMS Value= 0.180422 Quantization Noise Power= 0.032552');
disp('6 bits QEO: RMS Value= 0.045105 Quantization Noise Power= 0.002035');
disp('8 bits QEO: RMS Value= 0.011276 Quantization Noise Power= 0.000127');
disp('10 bits QEO: RMS Value= 0.002819 Quantization Noise Power= 0.000008');

x=[4 6 8 10];
y1=[0.180422 0.045105 0.011276 0.002819];
y2=[0.032552 0.002035 0.000127 0.000008];
plot(x,y1,'bo-');
hold on;
plot(x,y2,'ro-');
legend('RMS Value', 'Quantization Noise Power');
ylabel('Error');
xlabel('N-Bit Quantization');
title('Error vs N-Bit Quantization');

%% 2. Sin 10V
clc;
close all;
clear all;

t=0:1e-9:(2e-3)-(1e-9);
x=9.9999999*sin(2*pi*1e3*t);
quantile_interval=20/(2^4);
quantizer=floor(x/quantile_interval)*quantile_interval+(quantile_interval/2);

figure;
subplot(2,1,1);
plot(t,x);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Orignal Signal');
grid;
axis([0 2e-3 -12 12]);

subplot(2,1,2);
plot(t,quantizer);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Quantized Signal');
grid;
axis([0 2e-3 -12 12]);

figure;
error=x-quantizer;
plot(t,error);
xlabel('Time(s)');
ylabel('Error Voltage(V)');
title('Error');
axis([0 2e-3 -0.75 0.75]);
rms=sqrt(mean(error.^2));

quantization_noise_power=(quantile_interval)^2/12;

fprintf('RMS Value= %f\nQuantization Noise Power= %f\n',rms, quantization_noise_power);

figure;
histogram(error);
xlabel('Error(V)');
ylabel('Occurences');
title('Histogram of Errors');

%% 2. Sin 5 V
clc;
close all;
clear all;

t=0:1e-9:(2e-3)-(1e-9);
x=4.99999999*sin(2*pi*1e3*t);
quantile_interval=10/(2^4);
quantizer=floor(x/quantile_interval)*quantile_interval+(quantile_interval/2);

figure;
subplot(2,1,1);
plot(t,x);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Orignal Signal');
grid;
axis([0 2e-3 -6 6]);

subplot(2,1,2);
plot(t,quantizer);
xlabel('Time(s)');
ylabel('Voltage(V)');
title('Quantized Signal');
grid;
axis([0 2e-3 -6 6]);

figure;
error=x-quantizer;
plot(t,error);
xlabel('Time(s)');
ylabel('Error Voltage(V)');
title('Error');
axis([0 2e-3 -0.75 0.75]);
rms=sqrt(mean(error.^2));

quantization_noise_power=(quantile_interval)^2/12;

fprintf('RMS Value= %f\nQuantization Noise Power= %f\n',rms, quantization_noise_power);


figure;
histogram(error);
xlabel('Error(V)');
ylabel('Occurences');
title('Histogram of Errors');
fprintf('4 bits QEO: RMS Value\n');