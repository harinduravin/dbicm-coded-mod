EbNoVec = (5.5:0.04:5.94)';      % Eb/No values (dB)
berEst = zeros(size(EbNoVec));

% Delay profile
delay_profile = [0 0 1 0 1];
seed = 124;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = dvbs2ldpc(2/3);
ldpcEncoder = comm.LDPCEncoder(p);
ldpcDecoder = comm.LDPCDecoder(p,'DecisionMethod','Soft decision','OutputValue','Whole codeword');
ldpcDecoderHard = comm.LDPCDecoder(p);

msg0 = logical(randi([0 1],size(p,2)-size(p,1),1)); % older message 
msg1 = logical(randi([0 1],size(p,2)-size(p,1),1));
msg2 = logical(randi([0 1],size(p,2)-size(p,1),1)); 
msg3 = logical(randi([0 1],size(p,2)-size(p,1),1)); % newest message

% LDPC Encoding and interleaving
encData0 = ldpcEncoder(msg0);
inter0 = randintrlv(int8(encData0),seed); % Interleave.

% LDPC Encoding and interleaving
encData1 = ldpcEncoder(msg1);
inter1 = randintrlv(int8(encData1),seed); % Interleave.

% LDPC Encoding and interleaving
encData2 = ldpcEncoder(msg2);
inter2 = randintrlv(int8(encData2),seed); % Interleave.

% LDPC Encoding and interleaving
encData3 = ldpcEncoder(msg3);
inter3 = randintrlv(int8(encData3),seed); % Interleave.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input to the Delay module
s_2_p0 = reshape(inter0, length(inter0)/5, 5).';
s_2_p1 = reshape(inter1, length(inter1)/5, 5).';
s_2_p2 = reshape(inter2, length(inter2)/5, 5).';
s_2_p3 = reshape(inter3, length(inter3)/5, 5).';

for n = 1:length(EbNoVec)
   
    % Reset initialization values
    numErrs = 0;
    numBits = 0;
    iteration = 0;
    parallel_data_old = zeros(size(s_2_p0));
    s_2_pllr = zeros(size(s_2_p0)) + Inf;
 
    while numErrs < 30000 && numBits < 1e7

        iteration = iteration +1;

        switch mod(iteration,4)
            case 0
                parallel_data = s_2_p0;
            case 1
                parallel_data = s_2_p1;
            case 2
                parallel_data = s_2_p2;
            otherwise
                parallel_data = s_2_p3;
        end

        % Performing delay according to delay profile
        mixed_matrix = parallel_data;
        mixed_matrix(delay_profile == 1, :) = parallel_data_old(delay_profile == 1, :);

        parallel_data_old = parallel_data;

        % Ready to modulate
        mixed_matrix_flat = mixed_matrix(:);

        % Mapper
        modulated = APSK_32_mapper(int8(mixed_matrix_flat));
        % modulated = dvbsapskmod(int8(mixed_matrix_flat),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);

        % Pass through AWGN channel
        rx_symbols = add_awgn_32apsk(modulated, EbNoVec(n),2/3,5);

        %%%%%%%%%% Receiver %%%%%%%%%%

        % Demodulation
        demodulated = APSK_32_demapper_optimized(rx_symbols, EbNoVec(n));
        % demodulated =dvbsapskdemod(rx_symbols,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/((2/3)*2*10^((EbNoVec(n)+10*log10(5))/10)),'UnitAveragePower',true);

        % Converting to serial to parallel
        s_2_p_rx = reshape(demodulated, 5, length(demodulated)/5);

        if iteration ~= 1
            % Reverse delay module
            undecoded = s_2_p_rx_old;
            undecoded(delay_profile == 1, :) = s_2_p_rx(delay_profile == 1, :);
    
            % Ready to decode
            undecoded_T = undecoded.';
            undecoded_flat = undecoded_T(:);
    
            % Deinterleaving
            undecoded_deintered = randdeintrlv(undecoded_flat,seed); % Deinterleave.
    
            % Decoding
            decodedllr = ldpcDecoder(undecoded_deintered);
            % rxBits = ldpcDecoderHard(undecoded_deintered);
            rxBits = decodedllr(1:43200)<0;
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Storing soft LLR values for feedback
            
            % Interleaving
            interllr = randintrlv(decodedllr,seed); % Interleave.
            
            % Serial to parallel
            s_2_pllr = reshape(interllr, length(interllr)/5, 5).';
            s_2_pllr(delay_profile == 0, :) = 0;


            % Implementing hard feedback

            % s_2_pllr(s_2_pllr>0) = Inf;
            % s_2_pllr(s_2_pllr<0) = -Inf;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            MSG=rxBits;
                
            % Calculate the number of bit errors
            switch mod(iteration,4)
                case 1
                    nErrors = biterr(msg0,MSG);
                case 2
                    nErrors = biterr(msg1,MSG);
                case 3
                    nErrors = biterr(msg2,MSG);
                otherwise
                    nErrors = biterr(msg3,MSG);
            end

            % Increment the error and bit counters
            numErrs = numErrs + nErrors;
            numBits = numBits + length(msg0);
        end

        demodulated_acc = APSK_32_feedback_demapper(rx_symbols, EbNoVec(n), s_2_pllr, delay_profile);

        % Converting to serial to parallel
        s_2_p_rx_old = reshape(demodulated_acc, 5, length(demodulated_acc)/5);

    end

    numErrs
    numBits

    % Estimate the BER
    berEst(n) = numErrs/numBits
end

semilogy(EbNoVec,berEst,'-.')
grid
legend('Estimated BER ldpc')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')










% test_sym = [0.964523701047841 + 0.675167092165942j
% 0.758301823069097 + 0.805935836425300i;
% -0.384151790899565 + 1.19932707834310i;
% -0.280453462114944 + 0.840049681092105i];
% % test_sym = [0.46762566+0.96762566j];
% 
% demodulated_0_m = APSK_32_demapper(test_sym, 6);
% demodulated_0 = dvbsapskdemod(test_sym,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/(2*10^((6+10*log10(5))/10)),'UnitAveragePower',true); 
% 
% 
% modulated_0 = dvbsapskmod(int8([0;1;0;1;1]),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);



% % Shannon limit for the scheme 4.3517 dB
% 
% % EbNoVec = (12.4:0.04:12.8)';      % Eb/No values (dB)
% % EbNoVec = (10.4:0.04:10.8)';      % Eb/No values (dB)
% % EbNoVec = (8.6:0.04:9)';      % Eb/No values (dB)
% EbNoVec = (5.0:0.04:6.44)';      % Eb/No values (dB)
% berEst = zeros(size(EbNoVec));
% % bpskModulator = comm.BPSKModulator;
% % bpskDEModulator = comm.BPSKDemodulator('DecisionMethod',"Log-likelihood ratio");
% p = dvbs2ldpc(2/3);
% ldpcEncoder = comm.LDPCEncoder(p);
% ldpcDecoder = comm.LDPCDecoder(p,'DecisionMethod','Soft decision','OutputValue','Whole codeword');
% ldpcDecoderHard = comm.LDPCDecoder(p);
% Ebno = 6.1;
% 
% msg0 = logical(randi([0 1],size(p,2)-size(p,1),1)); % older message 
% msg1 = logical(randi([0 1],size(p,2)-size(p,1),1));
% msg2 = logical(randi([0 1],size(p,2)-size(p,1),1)); 
% msg3 = logical(randi([0 1],size(p,2)-size(p,1),1)); % newest message
% 
% % LDPC Encoding and interleaving
% encData0 = ldpcEncoder(msg0);
% inter0 = randintrlv(int8(encData0),25689); % Interleave.
% 
% % LDPC Encoding and interleaving
% encData1 = ldpcEncoder(msg1);
% inter1 = randintrlv(int8(encData1),25689); % Interleave.
% 
% % LDPC Encoding and interleaving
% encData2 = ldpcEncoder(msg2);
% inter2 = randintrlv(int8(encData2),25689); % Interleave.
% 
% % LDPC Encoding and interleaving
% encData3 = ldpcEncoder(msg3);
% inter3 = randintrlv(int8(encData3),25689); % Interleave.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Delay module is added here
% 
% % Input to the Delay module
% s_2_p0 = reshape(inter0, length(inter0)/5, 5).';
% s_2_p1 = reshape(inter1, length(inter1)/5, 5).';
% s_2_p2 = reshape(inter2, length(inter2)/5, 5).';
% s_2_p3 = reshape(inter3, length(inter3)/5, 5).';
% 
% % Delay profile
% delay_profile = [1 0 1 0 1];
% 
% % Create the first mixed matrix
% mixed_matrix_0 = s_2_p0;
% mixed_matrix_0(delay_profile == 1, :) = 0;
% 
% % Create the older mixed matrix
% mixed_matrix_1 = s_2_p1;
% mixed_matrix_1(delay_profile == 1, :) = s_2_p0(delay_profile == 1, :);
% 
% % Create the middle mixed matrix
% mixed_matrix_2 = s_2_p2;
% mixed_matrix_2(delay_profile == 1, :) = s_2_p1(delay_profile == 1, :);
% 
% % Create the newer mixed matrix
% mixed_matrix_3 = s_2_p3;
% mixed_matrix_3(delay_profile == 1, :) = s_2_p2(delay_profile == 1, :);
% 
% % Ready to modulate
% mixed_matrix_flat_0 = mixed_matrix_0(:);
% mixed_matrix_flat_1 = mixed_matrix_1(:);
% mixed_matrix_flat_2 = mixed_matrix_2(:);
% mixed_matrix_flat_3 = mixed_matrix_3(:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modulated_0 = dvbsapskmod(int8(mixed_matrix_flat_0),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);
% modulated_1 = dvbsapskmod(int8(mixed_matrix_flat_1),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);
% modulated_2 = dvbsapskmod(int8(mixed_matrix_flat_2),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);
% modulated_3 = dvbsapskmod(int8(mixed_matrix_flat_3),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);
% 
% rx_symbols_0 = add_awgn_32apsk(modulated_0, Ebno,2/3,5);
% rx_symbols_1 = add_awgn_32apsk(modulated_1, Ebno,2/3,5);
% rx_symbols_2 = add_awgn_32apsk(modulated_2, Ebno,2/3,5);
% rx_symbols_3 = add_awgn_32apsk(modulated_3, Ebno,2/3,5);
% 
% demodulated_0 = APSK_32_demapper_optimized(rx_symbols_0, Ebno);
% demodulated_1 = APSK_32_demapper_optimized(rx_symbols_1, Ebno);
% demodulated_2 = APSK_32_demapper_optimized(rx_symbols_2, Ebno);
% demodulated_3 = APSK_32_demapper_optimized(rx_symbols_3, Ebno);
% 
% % demodulated_0 = dvbsapskdemod(rx_symbols_0,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/(2*10^((Ebno+10*log10(5))/10)),'UnitAveragePower',true); 
% % demodulated_1 = dvbsapskdemod(rx_symbols_1,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/(2*10^((Ebno+10*log10(5))/10)),'UnitAveragePower',true); 
% % demodulated_2 = dvbsapskdemod(rx_symbols_2,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/(2*10^((Ebno+10*log10(5))/10)),'UnitAveragePower',true); 
% % demodulated_3 = dvbsapskdemod(rx_symbols_3,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/(2*10^((Ebno+10*log10(5))/10)),'UnitAveragePower',true); 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inverse delay module is added here
% 
% s_2_p0_delayed = reshape(demodulated_0, 5, length(demodulated_0)/5);
% s_2_p1_delayed = reshape(demodulated_1, 5, length(demodulated_1)/5);
% s_2_p2_delayed = reshape(demodulated_2, 5, length(demodulated_2)/5);
% s_2_p3_delayed = reshape(demodulated_3, 5, length(demodulated_3)/5);
% 
% undecoded_0 = s_2_p0_delayed;
% undecoded_0(delay_profile == 1, :) = s_2_p1_delayed(delay_profile == 1, :);
% 
% undecoded_1 = s_2_p1_delayed;
% undecoded_1(delay_profile == 1, :) = s_2_p2_delayed(delay_profile == 1, :);
% 
% undecoded_2 = s_2_p2_delayed;
% undecoded_2(delay_profile == 1, :) = s_2_p3_delayed(delay_profile == 1, :);
% 
% % Ready to decode
% undecoded_0_T = undecoded_0.';
% undecoded_flat_0 = undecoded_0_T(:);
% 
% undecoded_1_T = undecoded_1.';
% undecoded_flat_1 = undecoded_1_T(:);
% 
% undecoded_2_T = undecoded_2.';
% undecoded_flat_2 = undecoded_2_T(:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% undecoded_deintered_0 = randdeintrlv(undecoded_flat_0,25689); % Deinterleave.
% undecoded_deintered_1 = randdeintrlv(undecoded_flat_1,25689); % Deinterleave.
% undecoded_deintered_2 = randdeintrlv(undecoded_flat_2,25689); % Deinterleave.
% 
% decodedllr0 = ldpcDecoder(undecoded_deintered_0);
% rxBits = ldpcDecoderHard(undecoded_deintered_1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Implementing the feedback demapper
% 
% % Interleaving
% interllr0 = randintrlv(decodedllr0,25689); % Interleave.
% 
% % Serial to parallel
% s_2_pllr0 = reshape(interllr0, length(interllr0)/5, 5).';
% 
% demodulated_1_improved = APSK_32_feedback_demapper(rx_symbols_1, Ebno, s_2_pllr0, delay_profile);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inverse delay module is added here
% s_2_p1_delayed_improved = reshape(demodulated_1_improved, 5, length(demodulated_1)/5);
% 
% undecoded_1_improved = s_2_p1_delayed_improved;
% undecoded_1_improved(delay_profile == 1, :) = s_2_p2_delayed(delay_profile == 1, :);
% 
% % Ready to decode
% undecoded_1_improved_T = undecoded_1_improved.';
% undecoded_improved_flat_1 = undecoded_1_improved_T(:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% undecoded_deintered_improved_1 = randdeintrlv(undecoded_improved_flat_1,25689); % Deinterleave.
% 
% rxBits_improved = ldpcDecoderHard(undecoded_deintered_improved_1);
% 
% chk = biterr(msg1,rxBits);
% chk_improved = biterr(msg1,rxBits_improved);
% 
% % Calculate the number of bit errors
% nErrors = biterr(msg2,MSG);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% for n = 1:length(EbNoVec)
% 
%     % Reset the error and bit counters
%     numErrs = 0;
%     numBits = 0;
% 
%     while numErrs < 30000 && numBits < 1e8
% 
%         % Mapper
%         dataMod_psk = APSK_32_mapper(inter);
%         dataMod_pskx = dvbsapskmod(int8(inter),32,'s2','2/3','InputType','bit','UnitAveragePower',true);
% 
%         % Pass through AWGN channel
%         dataMod_psk2 = add_awgn_32apsk(dataMod_psk, EbNoVec(n),2/3,5);
% 
%         % Demodulation
%         dataDeMod_psk3 = dvbsapskdemod(dataMod_psk2,32,'s2','2/3','OutputType','llr','NoiseVariance',1/(2*10^((EbNoVec(n)+10*log10(5))/10)));  
% 
%         % Deinterleaving
%         deinter = randdeintrlv(dataDeMod_psk3,25689); % Deinterleave.
% 
%         % Decoding
%         rxBits = ldpcDecoder(deinter);
% 
%         MSG=rxBits;
%         chk = isequal(msg,MSG);
% 
%         % Calculate the number of bit errors
%         nErrors = biterr(msg,MSG);
% 
%         % Increment the error and bit counters
%         numErrs = numErrs + nErrors;
%         numBits = numBits + length(msg);
%     end
% 
%     numErrs
%     numBits
% 
%     % Estimate the BER
%     berEst(n) = numErrs/numBits
% end
% 
% berTheory = berawgn(EbNoVec,'psk',2,'nondiff');
% semilogy(EbNoVec,berEst,'-.')
% grid
% legend('Estimated BER ldpc')
% xlabel('Eb/No (dB)')
% ylabel('Bit Error Rate')



function out_data = add_awgn_32apsk(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end