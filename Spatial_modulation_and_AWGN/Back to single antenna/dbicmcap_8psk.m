clear;
EbNoVec = (-15:2.5:20)';      % Eb/No values (dB)

N_bicm = 32000;
N = N_bicm*1/8;

% bits per symbol
m = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayscheme = zeros(1,m);

delayscheme(1) = 1;
delayscheme(3) = 1;

numDelaybits = sum(delayscheme(:) == 1);

delaylabels = cell(1, 2^numDelaybits);

for i = 0:(2^numDelaybits - 1)
    binary_seq = dec2bin(i, numDelaybits);
    delaylabels{i + 1} = binary_seq;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symbols = [ 
    1.00000000e+00+0.00000000e+00j,  7.07106781e-01+7.07106781e-01j,...
    6.12323400e-17+1.00000000e+00j, -7.07106781e-01+7.07106781e-01j,...
   -1.00000000e+00+1.22464680e-16j, -7.07106781e-01-7.07106781e-01j,...
   -1.83697020e-16-1.00000000e+00j,  7.07106781e-01-7.07106781e-01j
  ];

labels = {'000', '001', '011', '010', '110', '111', '101','100'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

msg_bicm = logical(randi([0 1],m,N_bicm));
symbmsgbicm = zeros(1,N_bicm);
for q = 1:N_bicm
    index = find(strcmp(labels, sprintf('%d',msg_bicm(:,q)')));
    symbmsgbicm(q) = symbols(index);
end

capEst = zeros(m,length(EbNoVec));
shannoncap = zeros(size(EbNoVec)); 
ccCap = zeros(size(EbNoVec)); 
bicmCap = zeros(m,length(EbNoVec)); 

for n = 1:length(EbNoVec)

    cccapval = 0;

    esno = 10^((EbNoVec(n)+10*log10(m))/10);
    noise_var = 1/(2*10^((EbNoVec(n)+10*log10(m))/10));
    shannoncap(n) = log2(1+esno);
 
    for j = 1:N_bicm/2

        modulated = pskmod(int8(msg_bicm(:,j)),8,0,'gray','InputType','bit');

        % Pass through AWGN channel
        y = add_awgn(modulated, EbNoVec(n),1,m);
        cccapval = cccapval + (1/(N_bicm/2))*(m-log2(prob(msg_bicm(1,j),0,y,EbNoVec(n),symbols,labels)/exp(eucl_dist_norm(modulated, y, noise_var))));
    end

    % Store the capacity
    ccCap(n) = cccapval;
end

for ii = 1:m
    for n = 1:length(EbNoVec)
        EbNoVec(n)
    
        llccap = 0;

        for j = 1:N_bicm
            
            % AWGN noise at receiver
            y = add_awgn(symbmsgbicm(j), EbNoVec(n),1,m);
    
            % Bit levels 
            bitlevels = msg_bicm(:,j);

            llccap = llccap + (1/N_bicm)*(log2(prob(bitlevels(ii),0,y,EbNoVec(n),symbols,labels)/prob(bitlevels(ii),ii,y,EbNoVec(n),symbols,labels)));
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
                    y = add_awgn(symbmsg(j,zz), EbNoVec(n),1,m);
        
                    % Bit levels
                    bitlevels = msg(:,j);
        
                    cap = cap + (1/N)*(log2(prob(bitlevels(ii),0,y,EbNoVec(n),delaysymbollist,delaysymlabels)/prob(bitlevels(ii),ii,y,EbNoVec(n),delaysymbollist,delaysymlabels)));
                end
                llccap = llccap + (1/length(delaylabels))*cap;
            end
            % Store the capacity
            capEst(ii,n) = 1-llccap;
        end
    end
end
altstyles = {'bo--'; 'rs-.'; 'kd-'; 'g+--'; 'kv-'; 'b--'; 'r-'; 'k-'; 'g--'; 'g-';'b-'; 'mv-'};
Legend=cell(m,1);
for iter=1:m
    Legend{iter}=strcat('Bit -', num2str(iter));
end

figure
for z = 1:m
    plot(EbNoVec,capEst(z,:),altstyles{z})
    hold on;
end
grid
legend(Legend)
xlabel('Eb/No (dB)')
ylabel('Bit capacity')

figure
plot(EbNoVec,sum(capEst,1),altstyles{9})
hold on;
plot(EbNoVec,shannoncap,'-.')
plot(EbNoVec,ccCap,'-.')
plot(EbNoVec,sum(bicmCap,1),'-.')
grid
ylim([0 4])
legend('DBICM [1,0,1]','Shannon capacity','Const. Constr. capacity','BICM')
xlabel('Eb/No (dB)')
ylabel('Spectral Efficiency (bits/s/Hz)')

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

    noise_var = 1/(2*10^((ebno+10*log10(3))/10));
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