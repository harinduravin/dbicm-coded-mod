
EbNoVec = (3.8:0.1:5.2)';      % Eb/No values (dB)
berEst = zeros(size(EbNoVec));

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

msg = logical(randi([0 1],size(pcmatrix,2)-size(pcmatrix,1),1)); % older message

% Transmit and receive LDPC coded signal data
encData = ldpcEncoder(msg);
inter = randintrlv(int8(encData),2); % Interleave.

for n = 1:length(EbNoVec)
   
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;
 
    while numErrs < 4000 && numBits < 1e8
        
        % Mapper
        modulated = APSK_8_mapper(inter);

        % Pass through AWGN channel
        rx_symbols = add_awgn_8apsk(modulated, EbNoVec(n),3/4,3);

        % Demodulation
        demodulated = APSK_8_demapper_optimized(rx_symbols, EbNoVec(n));

        deinter = randdeintrlv(demodulated,2); % Deinterleave.

        % Decoding
        decodedllr = ldpcDecoder(deinter);
        rxBits = decodedllr(1:size(msg))<0;

        MSG=rxBits;
        chk = isequal(msg,MSG);
            
       % Calculate the number of bit errors
        nErrors = biterr(msg,MSG);
        
        % Increment the error and bit counters
        numErrs = numErrs + nErrors;
        numBits = numBits + length(msg);
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