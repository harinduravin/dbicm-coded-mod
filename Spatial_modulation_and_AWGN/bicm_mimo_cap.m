% Shannon limit for the scheme 4.3517 dB

EbNoVec = (-15:2:20)';      % Eb/No values (dB)
capEst = zeros(size(EbNoVec));
N = 30000;
ii = 2;
m = 3;

p = dvbs2ldpc(2/3);
% ldpcEncoder = comm.LDPCEncoder(p);
% ldpcDecoder = comm.LDPCDecoder(p);

% 2x2 mimo channel
Nr = 2;
Nt = 2;


msg = logical(randi([0 1],N,m));

for n = 1:length(EbNoVec)
    EbNoVec(n)

    llccap = 0;
   
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;

    % accumulate=zeros(2,2,N/2);

 
    for j = 1:(N/Nt)

        H = 1/sqrt(2)*(randn(Nr,Nt) + 1i*randn(Nr,Nt));
        
        % Bit levels considered
        allbits = int8(msg(Nt*(j-1)+1:Nt*j,:)');
        bitlevels = reshape(allbits, [], 1);

        % Transmit symbols for 2x2 mimo
        x = pskmod(bitlevels,8,0,'gray','InputType','bit');
        r = H*x;

        % AWGN noise at receiver
        y = add_awgn(r, EbNoVec(n),1,1);
        % 
        llccap = llccap + (1/(N/Nt))*(1-log2(prob(bitlevels(ii),0,y,EbNoVec(n),m,H)/prob(bitlevels(ii),ii,y,EbNoVec(n),m,H)));
                
    end

    % Store the capacity
    capEst(n) = llccap;
end

plot(EbNoVec,capEst,'-.')
grid
legend('Estimated BER ldpc')
xlabel('Eb/No (dB)')
ylabel('Bit capacity')

function out_data = add_awgn(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end

function sumprob = prob(b,ii,y,ebno,m,H)
    
    symbols = { 
        1.00000000e+00+0.00000000e+00j,  7.07106781e-01+7.07106781e-01j,...
        6.12323400e-17+1.00000000e+00j, -7.07106781e-01+7.07106781e-01j,...
       -1.00000000e+00+1.22464680e-16j, -7.07106781e-01-7.07106781e-01j,...
       -1.83697020e-16-1.00000000e+00j,  7.07106781e-01-7.07106781e-01j
      };

    labels = {'000', '001', '011', '010', '110', '111', '101','100'};

    noise_var = 1/(2*10^((ebno+10*log10(1))/10));
    sumprob = get_sum(y, symbols, labels, noise_var, ii,b,H,m);

end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
end

% function eucl_dist = eucl_dist_norm(a, r, noise_var,H)
%     real_diff = real(z1) - real(z2);
%     imag_diff = imag(z1) - imag(z2);
%     squared_diff = real_diff.^2 + imag_diff.^2;
%     eucl_dist = -squared_diff / (2*noise_var);
% 
%     % eucl_dist = -(r-H*a)'*(r-H*a)/(2*noise_var);
% end

function sumval = get_sum(r, symbols, labels, noise_var, position, value,H,m)
    
    
    if position == 0
        selected_symbols = cell2mat(symbols);
    else
        if value == 0
            valuestr = '0';
        else
            valuestr = '1';
        end
        symbol_indices = get_indices(mod(position-1,m)+1, valuestr, char(labels));
        selected_symbols = cell2mat(symbols(symbol_indices));
    end

    % Use ndgrid to generate combinations
    if ceil(position/m) <= 1
        [Array1, Array2] = ndgrid(selected_symbols,cell2mat(symbols));
    elseif ceil(position/m) == 2
        [Array1, Array2] = ndgrid(cell2mat(symbols), selected_symbols);
    end

    % Combine the elements into pairs
    sel_hyp_symbols = [Array1(:), Array2(:)];
    Ha = H*sel_hyp_symbols.';
    r_Ha = r -Ha;
    eucl_norms = diag(diag((r_Ha'*r_Ha)));
    sumval = trace(exp(-eucl_norms/(2*noise_var)));

    % eucl_norms = eucl_dist_norm(Ha, r, noise_var);
    % 
    % sumval = sum(exp(eucl_norms));
end

