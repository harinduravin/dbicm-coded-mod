clear;
% EbNoVec = (-15:2.5:20)';      % Eb/No values (dB)

N_bicm = 2^19;

% bits per symbol
m = 6;

M = 2^m;
x = (0:M-1)';
symbols = qammod(x,M,UnitAveragePower=true,PlotConstellation=false);
symbols = real(symbols) - 1i*imag(symbols);
symbols = symbols';

% gray labelled 64-QAM
labels = {'000000','000001','000010','000011','000100','000101','000110','000111','001000','001001','001010','001011','001100',...
    '001101','001110','001111','010000','010001','010010','010011','010100','010101','010110','010111','011000','011001','011010',...
    '011011','011100','011101','011110','011111','100000','100001','100010','100011','100100','100101','100110','100111','101000',...
    '101001','101010','101011','101100','101101','101110','101111','110000','110001','110010','110011','110100','110101','110110',...
    '110111','111000','111001','111010','111011','111100','111101','111110','111111'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

capEst = zeros(m,length(1:(2^m-2)));
bicmCap = zeros(m,length(1:(2^m-2))); 

EbNo = 2 - 10*log10(m);

msg_bicm = logical(randi([0 1],m,N_bicm));
symbmsgbicm = zeros(1,N_bicm);
for q = 1:N_bicm
    index = find(strcmp(labels, sprintf('%d',msg_bicm(:,q)')));
    symbmsgbicm(q) = symbols(index);
end

% shannoncap = zeros(size(EbNoVec)); 
% ccCap = zeros(size(EbNoVec)); 

cccapval = 0;

esno = 10^((EbNo+10*log10(m))/10);
noise_var = 1/(2*10^((EbNo+10*log10(m))/10));
shannoncap = log2(1+esno)

for j = 1:N_bicm

    % modulated = pskmod(int8(msg_bicm(:,j)),8,0,'gray','InputType','bit');
    modulated = qammod(int8(msg_bicm(:,j)),M,'InputType','bit',UnitAveragePower=true);
    modulated = real(modulated) - 1i*imag(modulated);

    % Pass through AWGN channel
    y = add_awgn(modulated, EbNo,1,m);
    cccapval = cccapval + (1/(N_bicm))*(m-log2(prob(msg_bicm(1,j),0,y,EbNo,symbols,labels)/exp(eucl_dist_norm(modulated, y, noise_var))));
end

% for n = 1:length(EbNoVec)
for n = 1:(2^m-2)
    n
    tic

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    binaryString = dec2bin(n, m); % m specifies the minimum number of bits
    delayscheme = str2num(binaryString')'; % Convert binary string to binary vector
    
    numDelaybits = sum(delayscheme(:) == 1);
    
    delaylabels = cell(1, 2^numDelaybits);
    
    for i = 0:(2^numDelaybits - 1)
        binary_seq = dec2bin(i, numDelaybits);
        delaylabels{i + 1} = binary_seq;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = N_bicm*m/(2^numDelaybits);

    msg = logical(randi([0 1],m,N));
    
    symbmsg = zeros(N,2^numDelaybits);
    for i = 1:2^numDelaybits
        temporary_binary_matrix = msg;
        temporary_binary_matrix(logical(delayscheme),:) = repmat(str2num(char(delaylabels(i))'),1,N);
        for q = 1:N
            index = find(strcmp(labels, sprintf('%d',temporary_binary_matrix(:,q)')));
            symbmsg(q,i) = symbols(index);
        end
    end

    for ii = 1:m
    
        llccap = 0;

        for j = 1:N_bicm
            
            % AWGN noise at receiver
            y = add_awgn(symbmsgbicm(j), EbNo,1,m);
    
            % Bit levels 
            bitlevels = msg_bicm(:,j);

            llccap = llccap + (1/N_bicm)*(log2(prob(bitlevels(ii),0,y,EbNo,symbols,labels)/prob(bitlevels(ii),ii,y,EbNo,symbols,labels)));
        end

        bicmCap(ii,n) = 1-llccap;
        if delayscheme(ii) == 1        
            capEst(ii,n) = 1-llccap;
        end

        llccap = 0;

        if delayscheme(ii) ~= 1
            for zz = 1:length(delaylabels)
                delayedbits = delaylabels(zz);
                indexlist = get_dbicmindices(delayscheme, char(delayedbits), char(labels));

                delaysymlabels = labels(indexlist);
                delaysymbollist = symbols(:,indexlist);

                cap = 0;

                for j = 1:(N)

                    % AWGN noise at receiver
                    y = add_awgn(symbmsg(j,zz), EbNo,1,m);
        
                    % Bit levels
                    bitlevels = msg(:,j);
        
                    cap = cap + (1/N)*(log2(prob(bitlevels(ii),0,y,EbNo,delaysymbollist,delaysymlabels)/prob(bitlevels(ii),ii,y,EbNo,delaysymbollist,delaysymlabels)));
                end
                llccap = llccap + (1/length(delaylabels))*cap;
            end
            % Store the capacity
            capEst(ii,n) = 1-llccap;
        end
    end
    toc
end

sum(capEst,1)
sum(bicmCap,1)

% altstyles = {'bo--'; 'rs-.'; 'kd-'; 'g+--'; 'kv-'; 'b--'; 'r-'; 'k-'; 'g--'; 'g-';'b-'; 'mv-'};
% Legend=cell(m,1);
% for iter=1:m
%     Legend{iter}=strcat('Bit -', num2str(iter));
% end

% figure
% for z = 1:m
%     plot(EbNoVec,capEst(z,:),altstyles{z})
%     hold on;
% end
% grid
% legend(Legend)
% xlabel('Eb/No (dB)')
% ylabel('Bit capacity')
% 
% figure
% plot(EbNoVec,sum(capEst,1),altstyles{9})
% hold on;
% plot(EbNoVec,shannoncap,'-.')
% plot(EbNoVec,ccCap,'-.')
% plot(EbNoVec,sum(bicmCap,1),'-.')
% grid
% ylim([0 4])
% legend('DBICM [1,0,1,0]','Shannon capacity','Const. Constr. capacity','BICM')
% xlabel('Eb/No (dB)')
% ylabel('Spectral Efficiency (bits/s/Hz)')

function out_data = add_awgn(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end

function sumprob = prob(b,ii,y,ebno,symbollist,symlabels)

    noise_var = 1/(2*10^((ebno+10*log10(6))/10));
    sumprob = get_sum(y, symbollist, symlabels, noise_var, ii,b);

end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
end

function indices = get_dbicmindices(delayscheme, delayedbits, labels)
    indices = find(all(labels(:, logical(delayscheme)) == delayedbits,2));
end

function sumval = get_sum(complex_number, symbols, labels, noise_var, position, value)

    if position == 0
        selected_symbols = symbols;
    else
        if value == 0
            valuestr = '0';
        else
            valuestr = '1';
        end
        symbol_indices = get_indices(position, valuestr, char(labels));
        selected_symbols = symbols(:,symbol_indices);
    end

    % Combine the elements into pairs
    eucl_norms = vecnorm(complex_number - selected_symbols,2,1).^2;
    sumval = sum(exp(-eucl_norms/(2*noise_var)));

end

function eucl_dist = eucl_dist_norm(z1, z2, noise_var)
    real_diff = real(z1) - real(z2);
    imag_diff = imag(z1) - imag(z2);
    squared_diff = real_diff.^2 + imag_diff.^2;
    eucl_dist = -squared_diff / (2*noise_var);
end