% Shannon limit for the scheme 4.3517 dB

EbNoVec = (-15:2.5:20)';      % Eb/No values (dB)
capEst = zeros(size(EbNoVec));
N = 30000;
ii = 2;

p = dvbs2ldpc(2/3);
% ldpcEncoder = comm.LDPCEncoder(p);
% ldpcDecoder = comm.LDPCDecoder(p);
msg = logical(randi([0 1],N,5));


for n = 1:length(EbNoVec)
    EbNoVec(n)

    llccap = 0;
   
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;

 
    for j = 1:N
        
        % Transmit and receive LDPC coded signal data
        % encData = ldpcEncoder(msg);
        % inter = randintrlv(int8(msg),25689); % Interleave.
     
        modulated = dvbsapskmod(msg(j,:)',32,'s2x','2/3','InputType','bit','UnitAveragePower',true);

        % Pass through AWGN channel
        y = add_awgn_32apsk(modulated, EbNoVec(n),2/3,5);

        llccap = llccap + (1/N)*(1-log2(prob(msg(j,ii),0,y,EbNoVec(n))/prob(msg(j,ii),ii,y,EbNoVec(n))));
        
    end

    % Store the capacity
    capEst(n) = llccap;
end

plot(EbNoVec,capEst,'-.')
grid
legend('Estimated BER ldpc')
xlabel('Eb/No (dB)')
ylabel('Bit capacity')

function out_data = add_awgn_32apsk(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end

function sumprob = prob(b,ii,y,ebno)
    
    symbols = { 
        0.164079178283360 + 0.164079178283360i,...
        -0.164079178283360 + 0.164079178283360i,...
        -0.164079178283360 - 0.164079178283360i,...
        0.164079178283360 - 0.164079178283360i,...
        0.638788528436365 + 0.171162870328789i,...
        0.467625658107576 + 0.467625658107576i,...
        0.171162870328789 + 0.638788528436365i,...
        -0.171162870328789 + 0.638788528436365i,...
        -0.467625658107576 + 0.467625658107576i,...
        -0.638788528436365 + 0.171162870328789i,...
        -0.638788528436365 - 0.171162870328789i,...
        -0.467625658107576 - 0.467625658107576i,...
        -0.171162870328790 - 0.638788528436365i,...
        0.171162870328789 - 0.638788528436365i,...
        0.467625658107576 - 0.467625658107576i,...
        0.638788528436365 - 0.171162870328790i,...
        1.26309318727039 + 0.251244856101071i,...
        1.07079869947673 + 0.715484816502116i,...
        0.715484816502116 + 1.07079869947673i,...
        0.251244856101071 + 1.26309318727039i,...
        -0.251244856101070 + 1.26309318727039i,...
        -0.715484816502115 + 1.07079869947673i,...
        -1.07079869947673 + 0.715484816502116i,...
        -1.26309318727039 + 0.251244856101071i,...
        -1.26309318727039 - 0.251244856101071i,...
        -1.07079869947673 - 0.715484816502115i,...
        -0.715484816502116 - 1.07079869947673i,...
        -0.251244856101071 - 1.26309318727039i,...
        0.251244856101071 - 1.26309318727039i,...
        0.715484816502116 - 1.07079869947673i,...
        1.07079869947673 - 0.715484816502116i,...
        1.26309318727039 - 0.251244856101071i
      };

    labels = {'01111', '01101', '11101', '11111', '01110', '00110', '00111',...
               '00101', '00100', '01100', '11100', '10100', '10101', '10111',...
               '10110', '11110', '01011', '01010', '00010', '00011', '00001',...
               '00000', '01000', '01001', '11001', '11000', '10000', '10001',...
               '10011', '10010', '11010', '11011'};

    noise_var = 1/(2*10^((ebno+10*log10(5))/10));
    sumprob = get_sum(y, symbols, labels, noise_var, ii,b);

end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
end

function eucl_dist = eucl_dist_norm(z1, z2, noise_var)
    real_diff = real(z1) - real(z2);
    imag_diff = imag(z1) - imag(z2);
    squared_diff = real_diff.^2 + imag_diff.^2;
    eucl_dist = -squared_diff / (2*noise_var);
end

function sumval = get_sum(complex_number, symbols, labels, noise_var, position, value)
    
    if position == 0
        selected_symbols = cell2mat(symbols);
    else
        if value == 0
            valuestr = '0';
        else
            valuestr = '1';
        end
        symbol_indices = get_indices(position, valuestr, char(labels));
        selected_symbols = cell2mat(symbols(symbol_indices));
    end
    eucl_norms = eucl_dist_norm(selected_symbols, complex_number, noise_var);

    sumval = sum(exp(eucl_norms));
end
% 
% function sum = get_sum(complex_number, symbols, labels, noise_var, position)
%     llr_1 = get_indices(position, '1', char(labels));
%     llr_0 = get_indices(position, '0', char(labels));
% 
%     symbols_llr_1 = cell2mat(symbols(llr_1));
%     symbols_llr_0 = cell2mat(symbols(llr_0));
% 
%     eucl_llr_1 = eucl_dist_norm(symbols_llr_1, complex_number, noise_var);
%     eucl_llr_0 = eucl_dist_norm(symbols_llr_0, complex_number, noise_var);
% 
%     % Max log approximation
%     % llr = max(eucl_llr_1) - max(eucl_llr_0);
% 
%     % Accurate calculation
%     sum = log(sum(exp(eucl_llr_1)));
% end

