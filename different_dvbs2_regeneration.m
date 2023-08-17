% Shannon limit for the scheme 4.3517 dB

% EbNoVec = (12.4:0.04:12.8)';      % Eb/No values (dB)
% EbNoVec = (10.4:0.04:10.8)';      % Eb/No values (dB)
% EbNoVec = (8.6:0.04:9)';      % Eb/No values (dB)
EbNoVec = (11.2:0.04:11.6)';      % Eb/No values (dB)
berEst = zeros(size(EbNoVec));
code_rate = 2/3;
bits = 1;
% bpskModulator = comm.BPSKModulator;
% bpskDEModulator = comm.BPSKDemodulator('DecisionMethod',"Log-likelihood ratio");
p = dvbs2ldpc(code_rate);
ldpcEncoder = comm.LDPCEncoder(p);
ldpcDecoder = comm.LDPCDecoder(p);
msg = logical(randi([0 1],size(p,2)-size(p,1),1));
for n = 1:length(EbNoVec)
   
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;
 
    while numErrs < 30000 && numBits < 1e8
        
         % Transmit and receive LDPC coded signal data
        encData = ldpcEncoder(msg);
        inter = randintrlv(int8(encData),25689); % Interleave.
     
        % dataMod_psk = APSK_32_mapper(int8(inter));
        dataMod_psk = dvbsapskmod(int8(inter),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);

        % Pass through AWGN channel
        dataMod_psk2 = add_awgn_32apsk(dataMod_psk, EbNoVec(n),2/3,bits);
        
        dataDeMod_psk3 =dvbsapskdemod(dataMod_psk2,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/((code_rate)*10^((EbNoVec(n)+10*log10(bits))/10)),'UnitAveragePower',true); 
        % dataDeMod_psk3 = APSK_32_demapper_optimized(dataMod_psk2, EbNoVec(n)); 

        deinter = randdeintrlv(dataDeMod_psk3,25689); % Deinterleave.

        rxBits = ldpcDecoder(deinter);

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

function out_data = add_awgn_32apsk(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end