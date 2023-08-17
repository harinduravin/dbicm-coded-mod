
modulated = ones(100000,1)+ ones(100000,1)*1i;
ebno = 5;

% Pass through AWGN channel
rx_symbols = add_awgn_32apsk(modulated, ebno,2/3,5);

variance = var(rx_symbols);

theoretical = 1/((2/3)*2*10^((ebno+10*log10(5))/10));

codeRate = 2/3;
UncodedEbNo = 5;
CodedEbNo = UncodedEbNo + 10*log10(codeRate);
channel = comm.AWGNChannel(BitsPerSymbol=5,EbNo=CodedEbNo);
rxSig = channel(modulated);

matlabvar = var(real(rxSig));

function out_data = add_awgn_32apsk(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end