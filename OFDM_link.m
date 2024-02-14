clear
close all
clc
% Multicarrier and Multiantenna course project
%% Start by definining the system parameters.
%Active subcarriers
number_of_active_subcarriers = 1200;
%Subcarrier Spacing 15 KHz
subcarrier_spacing = 15e3;
%FFT size 
FFT_size = power(2, ceil(log2(number_of_active_subcarriers))); % FFT-size is 2048
%Sampling Frequency
Fs = subcarrier_spacing * FFT_size;  % sampling frequency
QAM_order_data = 64; % data modulation order, chosen based on my calculations
QAM_order_pilot = 2; % modulation order for training symbols
nB_data = log2(QAM_order_data); % number of bits per data symbol
nB_pilot = log2(QAM_order_pilot); % number of bits per pilot symbol
number_of_frames = 200;
number_of_pilots_per_frame = 1;
number_of_data_per_frame = 15;
frame_length = number_of_pilots_per_frame + number_of_data_per_frame; % number of symbols per frame


%% Generate the data streams.
tx_data = randi([0, QAM_order_data-1], number_of_active_subcarriers, ...
    number_of_data_per_frame * number_of_frames);
tx_pilot = randi([0, QAM_order_pilot-1], number_of_active_subcarriers, ...
    number_of_pilots_per_frame * number_of_frames);

%% Transmitter Implementation
% Generate the QAM symbols.
tx_data_symbols = qammod(tx_data, QAM_order_data, 'UnitAveragePower', true);
tx_pilot_symbols = qammod(tx_pilot, QAM_order_pilot, 'UnitAveragePower', true);
% Create the frame structure. The for loop allows to easily extend beyond two frames.
tx_frame = zeros(number_of_active_subcarriers, frame_length * number_of_frames);
for i = 0:number_of_frames-1
    tx_frame(:,i * frame_length + 1:(i + 1) * frame_length) = [tx_pilot_symbols(:, ...
        i * number_of_pilots_per_frame + 1 : ...
        (i + 1) * number_of_pilots_per_frame) ...
        tx_data_symbols(:, i * number_of_data_per_frame + 1 : ...
        (i + 1) * number_of_data_per_frame)];
end

% Zero-pad the unused subcarriers and shift zero-frequency component to center of spectrum.
% zero-padded symbols
zp_length = (FFT_size - number_of_active_subcarriers);
tx_zp_frame = [zeros(zp_length / 2, frame_length * number_of_frames);
              tx_frame;
              zeros(zp_length / 2, frame_length * number_of_frames)];
tx_szp_frame = fftshift(tx_zp_frame, 1); % shifted and zero-padded symbols

% Generate the OFDM signal using IFFT.
tx_bb_signal = ifft(tx_szp_frame, FFT_size); % baseband signal

% Add the cyclic prefix.
cp_duration = round(12e-6*Fs);
cp = tx_bb_signal(end-cp_duration+1:end,:);
tx_cp_bb_signal = [cp; tx_bb_signal]; % cyclic-prefixed baseband signal

%% Channel & Receiver Implementation
% Transmit over a multipath channel, whereas the channel changes for every frame.
% Channel model: Veh-AE
fd = 0; % Doppler spread, zero in block-fading model
gains = [0 -1 -9 -10 -15 -20]; % in dB
delays = [0 310 710 1090 1730 10000]*1e-9;% in seconds

chan = comm.RayleighChannel('SampleRate',Fs, ...
    'PathDelays',delays, ...
    'AveragePathGains',gains, ...
    'MaximumDopplerShift',fd);


 
snrs_dB = 0:2:50; % SNR range in dB

j = 0;
for snr_dB = snrs_dB
     for i = 0:number_of_frames-1

         % Pass the current frame through the channel
         tx_cp_bb_frame_signal = tx_cp_bb_signal(:, i * frame_length + 1:(i + 1) * frame_length);
         rx_cp_bb_signal = chan(tx_cp_bb_frame_signal(:));
         rx_cp_bb_signal = reshape(rx_cp_bb_signal, size(tx_cp_bb_frame_signal));
         release(chan); % creating a new channel instance for next frame

         % Add the noise
         rx_cp_bb_signal = awgn(rx_cp_bb_signal, snr_dB, 'measured');
    
         % Remove the cyclix prefix
         rx_bb_signal = rx_cp_bb_signal(cp_duration+1:end,:);
    
         % Perform OFDM demodulation
         rx_szp_symbols = fft(rx_bb_signal, FFT_size);
    
         % Remove inactive subcarriers
         rx_zp_symbols = ifftshift(rx_szp_symbols, 1);
         rx_symbols = rx_zp_symbols(1 + zp_length / 2:end - zp_length / 2,:);
    
         % Estimate the channel using only the pilot of the frame
         zf_coefficients = rx_symbols(:, 1) ./ tx_pilot_symbols(:, i + 1);
    
         % Equalize all received symbols using the coefficients
         rx_eq_symbols = rx_symbols ./ zf_coefficients;
    
         % Collect the equalized data symbols for postprocessing later
         rx_eq_symbols_all(:, i * number_of_data_per_frame + 1 : ...
         (i + 1) * number_of_data_per_frame) = rx_eq_symbols(:, 2:end);    
     end
     
   %% Demodulate the symbols.
   rx_data = qamdemod(rx_eq_symbols_all, QAM_order_data, 'UnitAveragePower', true);

   % Calculate the bit error rate with biterr()
   [numErrors, ber] = biterr(tx_data, rx_data);

   % Save the bit error rate into a bit error rates vector
   bitErrorRates(j+1) = ber;

   j = j + 1;
end


semilogy(snrs_dB, bitErrorRates, 'r', 'LineWidth', 1);
grid on
xlabel('SNR');
ylabel('BER');
title('Simulation Results')
%end


