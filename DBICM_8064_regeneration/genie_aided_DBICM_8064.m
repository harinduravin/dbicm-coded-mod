EbNoVec = (3.8:0.1:5.1)';      % Eb/No values (dB)
berEst = zeros(size(EbNoVec));

% Delay profile
delay_profile = [1 0 0];
m = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
blockSize = 336;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ldpcEncoder = comm.LDPCEncoder(pcmatrix);
ldpcDecoder = comm.LDPCDecoder(pcmatrix,'DecisionMethod','Soft decision','OutputValue','Whole codeword');
ldpcDecoderHard = comm.LDPCDecoder(pcmatrix);

msg0 = logical(randi([0 1],size(pcmatrix,2)-size(pcmatrix,1),1)); % older message 
msg1 = logical(randi([0 1],size(pcmatrix,2)-size(pcmatrix,1),1));
msg2 = logical(randi([0 1],size(pcmatrix,2)-size(pcmatrix,1),1)); 
msg3 = logical(randi([0 1],size(pcmatrix,2)-size(pcmatrix,1),1)); % newest message

% LDPC Encoding and interleaving
encData0 = ldpcEncoder(msg0);
inter0 = randintrlv(int8(encData0),2); % Interleave.

% LDPC Encoding and interleaving
encData1 = ldpcEncoder(msg1);
inter1 = randintrlv(int8(encData1),2); % Interleave.

% LDPC Encoding and interleaving
encData2 = ldpcEncoder(msg2);
inter2 = randintrlv(int8(encData2),2); % Interleave.

% LDPC Encoding and interleaving
encData3 = ldpcEncoder(msg3);
inter3 = randintrlv(int8(encData3),2); % Interleave.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input to the Delay module
s_2_p0 = reshape(inter0, length(inter0)/m, m).';
s_2_p1 = reshape(inter1, length(inter1)/m, m).';
s_2_p2 = reshape(inter2, length(inter2)/m, m).';
s_2_p3 = reshape(inter3, length(inter3)/m, m).';

for n = 1:length(EbNoVec)
   
    % Reset initialization values
    numErrs = 0;
    numBits = 0;
    iteration = 0;
    parallel_data_old = zeros(size(s_2_p0));
    s_2_pllr = zeros(size(s_2_p0)) + Inf;
 
    while numErrs < 4000 && numBits < 1e8

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
        modulated = APSK_8_mapper(int8(mixed_matrix_flat));

        % Pass through AWGN channel
        rx_symbols = add_awgn_8apsk(modulated, EbNoVec(n),3/4,3);

        %%%%%%%%%% Receiver %%%%%%%%%%

        % Demodulation
        demodulated = APSK_8_demapper_optimized(rx_symbols, EbNoVec(n));

        % Converting to serial to parallel
        s_2_p_rx = reshape(demodulated, m, length(demodulated)/m);

        if iteration ~= 1
            % Reverse delay module
            undecoded = s_2_p_rx_old;
            undecoded(delay_profile == 1, :) = s_2_p_rx(delay_profile == 1, :);
    
            % Ready to decode
            undecoded_T = undecoded.';
            undecoded_flat = undecoded_T(:);
    
            % Deinterleaving
            undecoded_deintered = randdeintrlv(undecoded_flat,2); % Deinterleave.
    
            % Decoding
            decodedllr = ldpcDecoder(undecoded_deintered);
            % rxBits = ldpcDecoderHard(undecoded_deintered);
            rxBits = decodedllr(1:size(msg0))<0;
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Storing soft LLR values for feedback
            
            % Interleaving
            interllr = randintrlv(decodedllr,2); % Interleave.
            
            % Serial to parallel
            s_2_pllr = reshape(interllr, length(interllr)/m, m).';
            s_2_pllr(delay_profile == 0, :) = 0;


            % Implementing hard feedback

            s_2_pllr(s_2_pllr>0) = Inf;
            s_2_pllr(s_2_pllr<0) = -Inf;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            MSG=rxBits;
                
            % Calculate the number of bit errors
            switch mod(iteration,4)
                case 1
                    nErrors = biterr(msg0,MSG);
                    s_2_pllr(s_2_p0 == 0) = Inf;
                    s_2_pllr(s_2_p0 == 1) = -Inf;
                case 2
                    nErrors = biterr(msg1,MSG);
                    s_2_pllr(s_2_p1 == 0) = Inf;
                    s_2_pllr(s_2_p1 == 1) = -Inf;
                case 3
                    nErrors = biterr(msg2,MSG);
                    s_2_pllr(s_2_p2 == 0) = Inf;
                    s_2_pllr(s_2_p2 == 1) = -Inf;
                otherwise
                    nErrors = biterr(msg3,MSG);
                    s_2_pllr(s_2_p3 == 0) = Inf;
                    s_2_pllr(s_2_p3 == 1) = -Inf;
            end

            % Increment the error and bit counters
            numErrs = numErrs + nErrors;
            numBits = numBits + length(msg0);
        end

        demodulated_acc = APSK_8_feedback_demapper(rx_symbols, EbNoVec(n), s_2_pllr, delay_profile);

        % Converting to serial to parallel
        s_2_p_rx_old = reshape(demodulated_acc, m, length(demodulated_acc)/m);

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

function out_data = add_awgn_8apsk(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end