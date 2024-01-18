clear;
EbNoVec = (-15:2.5:20)';      % Eb/No values (dB)
N = 32000;

capEst = zeros(N,length(EbNoVec));
shannoncap = zeros(size(EbNoVec)); 
ccCap = zeros(N,length(EbNoVec)); 

m = 3;

msg = logical(randi([0 1],N,m));


for n = 1:length(EbNoVec)
    EbNoVec(n)

    llccap = zeros(N,1);
    cccapval = zeros(N,1);

    esno = 10^((EbNoVec(n)+10*log10(m))/10);
    noise_var = 1/(2*10^((EbNoVec(n)+10*log10(m))/10));
    shannoncap(n) = log2(1+esno);

    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;

 
    for j = 1:N
        
        % Transmit and receive LDPC coded signal data
        % encData = ldpcEncoder(msg);
        % inter = randintrlv(int8(msg),25689); % Interleave.
     
        modulated = pskmod(int8(msg(j,:)'),8,0,'gray','InputType','bit');

        % Pass through AWGN channel
        y = add_awgn_32apsk(modulated, EbNoVec(n),1,m);
        cccapval(j) = (m-log2(prob(msg(j,1),0,y,EbNoVec(n),m)/exp(eucl_dist_norm(modulated, y, noise_var))));

        llccap(j) = (1-log2(prob(msg(j,1),0,y,EbNoVec(n),m)/prob(msg(j,1),1,y,EbNoVec(n),m))) + (1-log2(prob(msg(j,2),0,y,EbNoVec(n),m)/prob(msg(j,2),2,y,EbNoVec(n),m))) + (1-log2(prob(msg(j,3),0,y,EbNoVec(n),m)/prob(msg(j,3),3,y,EbNoVec(n),m)));
    end

    % Store the capacity
    capEst(:,n) = llccap;

    ccCap(:,n) = cccapval;
end

figure
plot(EbNoVec,mean(capEst),'-.')
hold on;
plot(EbNoVec,shannoncap,'-.')
plot(EbNoVec,mean(ccCap),'-.')
grid
legend('BICM capacity','Shannon capacity','Const. Constr. capacity')
xlabel('Eb/No (dB)')
ylim([0 4])
ylabel('Spectral Efficiency (bits/s/Hz)')

figure
plot(EbNoVec,std(capEst)/sqrt(N),'-.')
hold on;
plot(EbNoVec,std(ccCap)/sqrt(N),'-.')
grid
legend('BICM capacity std','Const. Constr. capacity std')
xlabel('Eb/No (dB)')
ylabel('Spectral Efficiency (bits/s/Hz)')

function out_data = add_awgn_32apsk(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end

function sumprob = prob(b,ii,y,ebno,m)
    
    symbols = { 
        1.00000000e+00+0.00000000e+00j,  7.07106781e-01+7.07106781e-01j,...
        6.12323400e-17+1.00000000e+00j, -7.07106781e-01+7.07106781e-01j,...
       -1.00000000e+00+1.22464680e-16j, -7.07106781e-01-7.07106781e-01j,...
       -1.83697020e-16-1.00000000e+00j,  7.07106781e-01-7.07106781e-01j
      };

    labels = {'000', '001', '011', '010', '110', '111', '101','100'};

    noise_var = 1/(2*10^((ebno+10*log10(m))/10));
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

