
%% Idea of the experiment
% The modulation of the binary stream can be done directly with help of trigonometric functions in case of QPSK with phase rotation $\pi/4$  and Gray mapping.
% We've supposed that this way will be faster than implemented in MatLab functions and objects.

%% Algorithm
% To implement this hint:
%
% # Split the bit stream into even and odd bits.
% # Apply to these groups NRZ (Non Return-to-Zero) procedure (e.g., multiply with -2 and add 1).
% # Multiply odd bits with j (make it an imaginary part of the envelope)
% # Multiply the groupes with $cos(\pi/4)$ .
% 
% To implement soft demodulation:
%
% # For odd bits select the real part of the received signal, and imaginary part for even bits - LLR calculation feature.
% # Apply the following formula $(x+cos(\pi / 4))^2 - (x-cos(\pi / 4))^2$ , where $x$ is the eithe real or imaginary part of the received signal.
% # Collect demodulated message.

%% Verification
% To verify the suggestion the following MatLab functions and objects were used:
%
%%
% * <https://www.mathworks.com/help/comm/ref/pskmod.html pskmod>
% * <https://www.mathworks.com/help/comm/ref/comm.qpskmodulator-system-object.html comm.QPSKModulator>
% * <https://www.mathworks.com/help/comm/ref/comm.qpskdemodulator-system-object.html?s_tid=doc_ta comm.QPSKDemodulator>
% 
% The message length: 10000 bits 
%
% Number of trials: 10000 
% 
% MatLab version: MatLab 2021a

clear all; clc; close all

hModulator = comm.QPSKModulator('BitInput', true,'PhaseOffset',pi/4);

hDemod = comm.QPSKDemodulator('PhaseOffset', pi/4, 'BitOutput', true,...
    'DecisionMethod','Log-likelihood ratio');
hDemod_2 = comm.QPSKDemodulator('PhaseOffset', pi/4, 'BitOutput', true,...
    'DecisionMethod', 'Approximate log-likelihood ratio'); 

s = randi(1,10000,1); % test message
d = zeros(length(s),1); 
s2 = step(hModulator,s); % modulation
noisy_qpsk = awgn(s2,10,'measured','dB'); %AWGN

fast_qpsk = @(s) ( (s(2:2:end)*(-2)+1) + j*(s(1:2:end-1)*(-2)+1) )*cos(pi/4);

d_2 = @(z) (imag(z) + cos(pi/4)).^2 - (imag(z) - cos(pi/4)).^2;
d_1 = @(z) (real(z) + cos(pi/4)).^2 - (real(z) - cos(pi/4)).^2;

for n = 1:10000

    tic;
    mm = fast_qpsk(s);
    t1(n) = toc;

    tic;
    a = step(hModulator,s);
    t2(n) = toc;

    tic;
    c = bi2de([s(1:2:end-1) s(2:2:end)],'left-msb');
    b = pskmod(c,4,pi/4,'gray');
    t3(n) = toc;

    tic;
    demod_msg = step(hDemod,noisy_qpsk);
    t4(n) = toc;

    tic;
    demod_msg_2 = step(hDemod_2,noisy_qpsk);
    t5(n) = toc;

    tic;
    d(1:2:end-1) = d_2(noisy_qpsk);
    d(2:2:end) = d_1(noisy_qpsk);
    t6(n) = toc;
end

%disp ('My algorithm (modulation)')
My_Algorithm_Modulation = mean(t1) % 8.8838e-05

%disp ('comm.QPSKModulator')
commQPSKModulator = mean(t2) % 1.4483e-04

%disp ('pskmod')
pskmod = mean(t3) % 6.3411e-04

%disp ('My algorithm (demodulation)')
My_Algorithm_Demodulation = mean(t6) % 7.2917e-05

%disp ('Exact LLR')
ExactLLR = mean(t4) % 0.0016

%disp ('Approximate LLR')
ApproximateLLR = mean(t5) % 4.4542e-04

% Do the algorithms provide the similar solutions?
% Yes. 
check1 = a - mm; % ->1e-16
check2 = demod_msg - demod_msg_2; % ->1e-16
check3 = demod_msg - d; % ->1e-16
