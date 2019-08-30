%% EECS 4214 Lab 5 Niruyan Rakulan 214343438
%% 10 db
% 1.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%given sequence
seq=[1 0 1 1 0 0 0 1];

%repeat 100 times
test=repmat(seq,1,100);

%map to voltages
v_test=test*10-5;

%geneate waveform
sample_time=10e-9;
t=0:sample_time:1e-6*length(test)-sample_time;
i=1e-6;
n=1;
for j=1:length(t)
    if t(j)<=n*i
        y(j)=v_test(n);
    else
        n=n+1;
        y(j)=v_test(n);
        
    end
end

plot(t,y);
title('Test Vector');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);


% 1.2 noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=10;
z=awgn(y,SNR,'measured');
figure;
subplot(3,1,1);
plot(t,y);
title('Original Signal');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);
subplot(3,1,2);
plot(t,z);
title(['Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);


% 1.3 Sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Since 100 samples per pulse, take 50 (i.e 50 , 150, 250 index of y) to get half
for j=1:length(test)
    s_time(j)=t((2*j-1)*50);
    sample(j)=z((2*j-1)*50);
end

subplot(3,1,3);
stem(s_time,sample);
title(['Sampled Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);

% 1.4 Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% threshold at 0
yo=0;
for j=1:length(test)
    if sample(j)>yo
        detect(j)=1;
    else
        detect(j)=0;
    end
end

errors=0;
for j=1:length(detect)
    if (detect(j)~=test(j))
        errors=errors+1;
    end
end

figure;
h=histogram(z);
morebins(h);
title(['PDF with SNR:',num2str(SNR),' dB']);
xlabel('Time(s)');
ylabel('Occurrences');

BER=errors/length(test);

fprintf('BER of %f, and Threshold of %f, at SNR of %i dB \n',BER,yo,SNR);
%% 8 db
% 1.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%given sequence
seq=[1 0 1 1 0 0 0 1];

%repeat 100 times
test=repmat(seq,1,100);

%map to voltages
v_test=test*10-5;

%geneate waveform
sample_time=10e-9;
t=0:sample_time:1e-6*length(test)-sample_time;
i=1e-6;
n=1;
for j=1:length(t)
    if t(j)<=n*i
        y(j)=v_test(n);
    else
        n=n+1;
        y(j)=v_test(n);
        
    end
end

plot(t,y);
title('Test Vector');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);


% 1.2 noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=8;
z=awgn(y,SNR,'measured');
figure;
subplot(3,1,1);
plot(t,y);
title('Original Signal');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);
subplot(3,1,2);
plot(t,z);
title(['Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);


% 1.3 Sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Since 100 samples per pulse, take 50 (i.e 50 , 150, 250 index of y) to get half
for j=1:length(test)
    s_time(j)=t((2*j-1)*50);
    sample(j)=z((2*j-1)*50);
end

subplot(3,1,3);
stem(s_time,sample);
title(['Sampled Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);

% 1.4 Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% threshold at 0
yo=0;
for j=1:length(test)
    if sample(j)>yo
        detect(j)=1;
    else
        detect(j)=0;
    end
end

errors=0;
for j=1:length(detect)
    if (detect(j)~=test(j))
        errors=errors+1;
    end
end

figure;
h=histogram(z);
morebins(h);
title(['PDF with SNR:',num2str(SNR),' dB']);
xlabel('Time(s)');
ylabel('Occurrences');

BER=errors/length(test);

fprintf('BER of %f, and Threshold of %f, at SNR of %i dB \n',BER,yo,SNR);

%% 6 db
% 1.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%given sequence
seq=[1 0 1 1 0 0 0 1];

%repeat 100 times
test=repmat(seq,1,100);

%map to voltages
v_test=test*10-5;

%geneate waveform
sample_time=10e-9;
t=0:sample_time:1e-6*length(test)-sample_time;
i=1e-6;
n=1;
for j=1:length(t)
    if t(j)<=n*i
        y(j)=v_test(n);
    else
        n=n+1;
        y(j)=v_test(n);
        
    end
end

plot(t,y);
title('Test Vector');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);


% 1.2 noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=6;
z=awgn(y,SNR,'measured');
figure;
subplot(3,1,1);
plot(t,y);
title('Original Signal');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);
subplot(3,1,2);
plot(t,z);
title(['Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -15 15]);


% 1.3 Sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Since 100 samples per pulse, take 50 (i.e 50 , 150, 250 index of y) to get half
for j=1:length(test)
    s_time(j)=t((2*j-1)*50);
    sample(j)=z((2*j-1)*50);
end

subplot(3,1,3);
stem(s_time,sample);
title(['Sampled Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -15 15]);

% 1.4 Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% threshold at 0
yo=0;
for j=1:length(test)
    if sample(j)>yo
        detect(j)=1;
    else
        detect(j)=0;
    end
end

errors=0;
for j=1:length(detect)
    if (detect(j)~=test(j))
        errors=errors+1;
    end
end

figure;
h=histogram(z);
morebins(h);
title(['PDF with SNR:',num2str(SNR),' dB']);
xlabel('Time(s)');
ylabel('Occurrences');

BER=errors/length(test);

fprintf('BER of %f, and Threshold of %f, at SNR of %i dB \n',BER,yo,SNR);

%% 4 db
% 1.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%given sequence
seq=[1 0 1 1 0 0 0 1];

%repeat 100 times
test=repmat(seq,1,100);

%map to voltages
v_test=test*10-5;

%geneate waveform
sample_time=10e-9;
t=0:sample_time:1e-6*length(test)-sample_time;
i=1e-6;
n=1;
for j=1:length(t)
    if t(j)<=n*i
        y(j)=v_test(n);
    else
        n=n+1;
        y(j)=v_test(n);
        
    end
end

plot(t,y);
title('Test Vector');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);


% 1.2 noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=4;
z=awgn(y,SNR,'measured');
figure;
subplot(3,1,1);
plot(t,y);
title('Original Signal');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);
subplot(3,1,2);
plot(t,z);
title(['Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -15 15]);


% 1.3 Sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Since 100 samples per pulse, take 50 (i.e 50 , 150, 250 index of y) to get half
for j=1:length(test)
    s_time(j)=t((2*j-1)*50);
    sample(j)=z((2*j-1)*50);
end

subplot(3,1,3);
stem(s_time,sample);
title(['Sampled Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -15 15]);

% 1.4 Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% threshold at 0
yo=0;
for j=1:length(test)
    if sample(j)>yo
        detect(j)=1;
    else
        detect(j)=0;
    end
end

errors=0;
for j=1:length(detect)
    if (detect(j)~=test(j))
        errors=errors+1;
    end
end

figure;
h=histogram(z);
morebins(h);
title(['PDF with SNR:',num2str(SNR),' dB']);
xlabel('Time(s)');
ylabel('Occurrences');

BER=errors/length(test);

fprintf('BER of %f, and Threshold of %f, at SNR of %i dB \n',BER,yo,SNR);

%% 2 db
% 1.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%given sequence
seq=[1 0 1 1 0 0 0 1];

%repeat 100 times
test=repmat(seq,1,100);

%map to voltages
v_test=test*10-5;

%geneate waveform
sample_time=10e-9;
t=0:sample_time:1e-6*length(test)-sample_time;
i=1e-6;
n=1;
for j=1:length(t)
    if t(j)<=n*i
        y(j)=v_test(n);
    else
        n=n+1;
        y(j)=v_test(n);
        
    end
end

plot(t,y);
title('Test Vector');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -10 10]);


% 1.2 noise%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=2;
z=awgn(y,SNR,'measured');
figure;
subplot(3,1,1);
plot(t,y);
title('Original Signal');
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -15 15]);
subplot(3,1,2);
plot(t,z);
title(['Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -15 15]);


% 1.3 Sample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Since 100 samples per pulse, take 50 (i.e 50 , 150, 250 index of y) to get half
for j=1:length(test)
    s_time(j)=t((2*j-1)*50);
    sample(j)=z((2*j-1)*50);
end

subplot(3,1,3);
stem(s_time,sample);
title(['Sampled Noisy Signal with SNR:',num2str(SNR),' dB'])
xlabel('Times(s)');
ylabel('Voltage(V)');
axis([min(t) max(t) -15 15]);

% 1.4 Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% threshold at 0
yo=0;
for j=1:length(test)
    if sample(j)>yo
        detect(j)=1;
    else
        detect(j)=0;
    end
end

errors=0;
for j=1:length(detect)
    if (detect(j)~=test(j))
        errors=errors+1;
    end
end

figure;
h=histogram(z);
morebins(h);
title(['PDF with SNR:',num2str(SNR),' dB']);
xlabel('Time(s)');
ylabel('Occurrences');

BER=errors/length(test);

fprintf('BER of %f, and Threshold of %f, at SNR of %i dB \n',BER,yo,SNR)

%% BER vs SNR

figure;
BER_A=[0.002500 0.005000 0.030000 0.062500 0.112500];
SNR_A=[10 8 6 4 2];
plot(SNR_A,log10(BER_A));
title('log(BER) vs SNR with Threshold at 0');
xlabel('SNR(dB)');
ylabel('log(BER)');