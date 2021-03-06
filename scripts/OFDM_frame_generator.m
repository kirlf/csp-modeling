clear all; close all; clc

M = 4; % e.g. QPSK 
N_inf = 16; % number of subcarriers (information symbols, actually) in the frame
fr_len = 32; % the length of our OFDM frame
N_pil = fr_len - N_inf - 5; % number of pilots in the frame
pilots = [1; j; -1; -j]; % pilots (BPSK, in fact)

nulls_idx = [1, 2, fr_len/2, fr_len-1, fr_len]; % indexes of nulls

idx_1_start = 4;
idx_1_end = fr_len/2 - 2;

idx_2_start = fr_len/2 + 2;
idx_2_end =  fr_len - 3;


inf_idx_1 = (floor(linspace(idx_1_start, idx_1_end, N_inf/2))).'; 
inf_idx_2 = (floor(linspace(idx_2_start, idx_2_end, N_inf/2))).';

inf_ind = [inf_idx_1; inf_idx_2]; % simple concatenation

inf_and_nulls_idx = union(inf_ind, nulls_idx); %concatenation and ascending sorting
pilot_idx = setdiff(1:fr_len, inf_and_nulls_idx); %numbers in range from 1 to frame length 
% that don't overlape with inf_and_nulls_idx vector

%% Pilots vector
% it should be very convinient to insert pilots if we prepare before "long-vector"
pilots_len_psudo = floor(N_pil/length(pilots)); % floor rounds value to lower integer
% - now we know how many full pilots vectors OFDM-frame consists

mat_1 = pilots*ones(1, pilots_len_psudo); % rank-one matrix - linear algebra trick
resh = reshape(mat_1, pilots_len_psudo*length(pilots),1); % vectorization - linear algebra trick

tail_len = fr_len  - N_inf - length(nulls_idx) ...
				- length(pilots)*pilots_len_psudo; 
tail = pilots(1:tail_len); % "tail" of pilots vector
vec_pilots = [resh; tail]; % completed pilots vector that frame consists

message = randi([0 M-1], N_inf, 1); % decimal information symbols

if M >= 16
	info_symbols = qammod(message, M, pi/4);
else
	info_symbols = pskmod(message, M, pi/4);
end 

%% Frame construction
frame = zeros(fr_len,1);
frame(pilot_idx) = vec_pilots;
frame(inf_ind) = info_symbols
