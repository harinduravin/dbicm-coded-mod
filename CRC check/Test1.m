EbNoVec = (5.9:0.04:6.94)';      % Eb/No values (dB)
berEst = zeros(size(EbNoVec));

crc8 = comm.CRCGenerator('Polynomial','z^8 + z^2 + z + 1', ...
    'InitialConditions',1,'DirectMethod',true,'FinalXOR',1);


% Delay profile
possible_delay_prof = [[1 0 0 0 0]; [0 1 0 0 0]; [0 0 1 0 0]; [0 0 0 1 0]; [0 0 0 0 1]];
delay_profile = possible_delay_prof(4,:);
seed = 124;

p = dvbs2ldpc(2/3);
ldpcEncoder = comm.LDPCEncoder(p);
ldpcDecoder = comm.LDPCDecoder(p,'DecisionMethod','Soft decision','OutputValue','Whole codeword');

msg0 = crc8(logical(randi([0 1],size(p,2)-size(p,1)-8,1))); % older message 
msg1 = crc8(logical(randi([0 1],size(p,2)-size(p,1)-8,1)));
msg2 = crc8(logical(randi([0 1],size(p,2)-size(p,1)-8,1))); 
msg3 = crc8(logical(randi([0 1],size(p,2)-size(p,1)-8,1))); % newest message

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

            for q = 1:5
                % Reverse delay module
                undecoded = s_2_p_rx_old;
                % undecoded(possible_delay_prof(4,:) == 1, :) = s_2_p_rx(possible_delay_prof(4,:) == 1, :);
                undecoded(possible_delay_prof(q,:) == 1, :) = s_2_p_rx(possible_delay_prof(q,:) == 1, :);
                % Ready to decode
                undecoded_T = undecoded.';
                undecoded_flat = undecoded_T(:);
    
                % Deinterleaving
                undecoded_deintered = randdeintrlv(undecoded_flat,seed); % Deinterleave.
    
                % Decoding
                decodedllr = ldpcDecoder(undecoded_deintered);
                % rxBits = ldpcDecoderHard(undecoded_deintered);
                rxBits = decodedllr(1:43200)<0;
    
                nErrors1 = biterr(msg0,rxBits);
    
                nErrors2 = biterr(msg1,rxBits);
    
                nErrors3 = biterr(msg2,rxBits);
    
                nErrors4 = biterr(msg3,rxBits);

                checksum_encode = crc8(rxBits(1:end-8));
                checksum = checksum_encode(end-8+1:end);
                expectedChecksum = rxBits(end-8+1:end);
                isequal(checksum,expectedChecksum)

            end
            
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

function out_data = add_awgn_32apsk(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end