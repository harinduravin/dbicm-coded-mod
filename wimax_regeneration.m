% Shannon limit for the scheme 4.3517 dB

% EbNoVec = (12.4:0.04:12.8)';      % Eb/No values (dB)
% EbNoVec = (10.4:0.04:10.8)';      % Eb/No values (dB)
% EbNoVec = (8.6:0.04:9)';      % Eb/No values (dB)
EbNoVec = (5:0.2:9)';      % Eb/No values (dB)
berEst = zeros(size(EbNoVec));
% bpskModulator = comm.BPSKModulator;
% bpskDEModulator = comm.BPSKDemodulator('DecisionMethod',"Log-likelihood ratio");
P = [3 0 -1 -1 2 0 -1 3 7 -1 1 1 -1 -1 -1 -1 1 0 -1 -1 -1 -1 -1 -1
-1 -1 1 -1 36 -1 -1 34 10 -1 -1 18 2 -1 3 0 -1 0 0 -1 -1 -1 -1 -1
-1 -1 12 2 -1 15 -1 40 -1 3 -1 15 -1 2 13 -1 -1 -1 0 0 -1 -1 -1 -1
-1 -1 19 24 -1 3 0 -1 6 -1 17 -1 -1 -1 8 39 -1 -1 -1 0 0 -1 -1 -1
20 -1 6 -1 -1 10 29 -1 -1 28 -1 14 -1 38 -1 -1 0 -1 -1 -1 0 0 -1 -1
-1 -1 10 -1 28 20 -1 -1 8 -1 36 -1 9 -1 21 45 -1 -1 -1 -1 -1 0 0 -1
35 25 -1 37 -1 21 -1 -1 5 -1 -1 0 -1 4 20 -1 -1 -1 -1 -1 -1 -1 0 0
-1 6 6 -1 -1 -1 4 -1 14 30 -1 3 36 -1 14 -1 1 -1 -1 -1 -1 -1 -1 0];

blockSize = 96;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);

cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);


msg = logical(randi([0 1],size(pcmatrix,2)-size(pcmatrix,1),1));
for n = 1:length(EbNoVec)
   
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;
 
    while numErrs < 30000 && numBits < 1e8
        
         % Transmit and receive LDPC coded signal data
        encData = ldpcEncode(msg,cfgLDPCEnc);
        inter = randintrlv(int8(encData),25689); % Interleave.

        % Add an additional bit for making it a multiple of 5
        inter(end+1) = 0;


     
         %channel
        % dataMod_psk = bpskModulator(encData);
        dataMod_psk = dvbsapskmod(int8(inter),32,'s2x','2/3','InputType','bit','UnitAveragePower',true);
        % Pass through AWGN channel
        % channel = comm.AWGNChannel("EbNo",EbNoVec(n),'BitsPerSymbol',5);
        % dataMod_psk2=channel(dataMod_psk);

        dataMod_psk2 = add_awgn_32apsk(dataMod_psk, EbNoVec(n),2/3,5);
        
        % dataDeMod_psk3=bpskDEModulator(dataMod_psk2);
        % dataDeMod_psk3=dvbsapskdemod(dataMod_psk2,32,'s2','3/4','OutputType','bit');
        dataDeMod_psk3 =dvbsapskdemod(dataMod_psk2,32,'s2x','2/3','OutputType','llr','NoiseVariance',1/(2*10^((EbNoVec(n)+10*log10(5))/10)),'UnitAveragePower',true);  


        dataDeMod_psk3(end) = [];

        deinter = randdeintrlv(dataDeMod_psk3,25689); % Deinterleave.

        % rxBits = ldpcDecode(single(-2*deinter+1),cfgLDPCDec,100);
        rxBits = ldpcDecode(deinter,cfgLDPCDec,100,"DecisionType","hard");

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
berTheory = berawgn(EbNoVec,'psk',2,'nondiff');
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