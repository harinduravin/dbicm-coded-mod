clear;
% EbNoVec = (-10:2.5:25)';     % Eb/No values (dB)
EsNoVec = (-2:1:16)'; 
N_bicm = 1000000;

% bits per symbol
m = 3;

EbNoVec = EsNoVec -10*log10(m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayscheme = [0,0,1];

numDelaybits = sum(delayscheme(:) == 1);

delaylabels = cell(1, 2^numDelaybits);

for i = 0:(2^numDelaybits - 1)
    binary_seq = dec2bin(i, numDelaybits);
    delaylabels{i + 1} = binary_seq;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 2^m;
x = (0:M-1)';
symbols = pammod(x,M,0,'gray')/sqrt(21);
% symbols = dvbsapskmod(x,M,'s2x','90/180',UnitAveragePower=true,PlotConstellation=false);
symbols = symbols';

% Gray labelled APSK-16
labels = {'000','001','010','011','100','101','110','111'};

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

msg_bicm = logical(randi([0 1],m,N_bicm));
symbmsgbicm = zeros(1,N_bicm);
for q = 1:N_bicm
    index = find(strcmp(labels, sprintf('%d',msg_bicm(:,q)')));
    symbmsgbicm(q) = symbols(index);
end

capEst = zeros(m,length(EbNoVec));
dbicm2 = zeros(m,length(EbNoVec)); 
dbicm3 = zeros(m,length(EbNoVec)); 
dbicm4 = zeros(m,length(EbNoVec)); 
shannoncap = zeros(size(EbNoVec)); 
ccCap = zeros(size(EbNoVec)); 
bicmCap = zeros(m,length(EbNoVec)); 

parfor n = 1:length(EbNoVec)

    cccapval = 0;

    esno = 10^((EbNoVec(n)+10*log10(m))/10);
    noise_var = 1/(10^((EbNoVec(n)+10*log10(m))/10));
    shannoncap(n) = log2(1+esno);
 
    for j = 1:N_bicm
        modulated = symbmsgbicm(j);
        % modulated = pskmod(int8(msg_bicm(:,j)),8,0,'gray','InputType','bit');
        % modulated = qammod(int8(msg_bicm(:,j)),16,'InputType','bit',UnitAveragePower=true);
        % modulated = dvbsapskmod(int8(msg_bicm(:,j)),M,'s2x','90/180',UnitAveragePower=true,PlotConstellation=false,InputType ='bit');

        % Pass through AWGN channel
        y = add_awgn(modulated, EbNoVec(n),1,m);
        cccapval = cccapval + (1/(N_bicm))*(m-log2(prob(msg_bicm(1,j),0,y,EbNoVec(n),symbols,labels)/exp(eucl_dist_norm(modulated, y, noise_var))));
    end

    % Store the capacity
    ccCap(n) = cccapval;
end

for ii = 1:m
    parfor n = 1:length(EbNoVec)
        EsNoVec(n)
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayscheme = [0,1,0];

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
    parfor n = 1:length(EbNoVec)
        EsNoVec(n)

        if delayscheme(ii) == 1        
            dbicm2(ii,n) = bicmCap(ii,n);
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
            dbicm2(ii,n) = 1-llccap;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayscheme = [1,0,1];

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
    parfor n = 1:length(EbNoVec)
        EsNoVec(n)

        if delayscheme(ii) == 1        
            dbicm3(ii,n) = bicmCap(ii,n);
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
            dbicm3(ii,n) = 1-llccap;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayscheme = [0,1,1];

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
    parfor n = 1:length(EbNoVec)
        EsNoVec(n)

        if delayscheme(ii) == 1        
            dbicm4(ii,n) = bicmCap(ii,n);
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
            dbicm4(ii,n) = 1-llccap;
        end
    end
end
altstyles = {'bo--'; 'rs-.'; 'kd-'; 'g+--'; 'kv-'; 'b--'; 'r-'; 'k-'; 'g--'; 'g-';'b-'; 'mv-'};
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
% xlabel('Es/No (dB)')
% ylabel('Bit capacity')

figure
plot(EsNoVec,2*sum(capEst,1),altstyles{9})
hold on;
plot(EsNoVec,2*sum(dbicm2,1),altstyles{8})
plot(EsNoVec,2*sum(dbicm3,1),altstyles{7})
plot(EsNoVec,2*sum(dbicm4,1),altstyles{6})
plot(EsNoVec,shannoncap,'-.')
plot(EsNoVec,2*ccCap,'-.')
plot(EsNoVec,2*sum(bicmCap,1),'-.')
grid
ylim([0 6.1])
legend('DBICM [0,0,1,0,0,1]','DBICM [0,1,0,0,1,0]','DBICM [1,0,1,1,0,1]','DBICM [0,1,1,0,1,1]','log_2(1+SNR)','CM capacity','BICM')
xlabel('Es/No (dB)')
ylabel('Spectral Efficiency (bits/s/Hz)')


speceff = linspace(0.8,4.9,10);
shannondb = 10*log10(2.^(speceff)-1);
bicmdb = interp1(2*sum(bicmCap,1),EsNoVec,speceff,"linear");
ccCapdb = interp1(2*ccCap,EsNoVec,speceff,"linear");
dbicmdb1 = interp1(2*sum(capEst,1),EsNoVec,speceff,"linear");
dbicmdb2 = interp1(2*sum(dbicm2,1),EsNoVec,speceff,"linear");
dbicmdb3 = interp1(2*sum(dbicm3,1),EsNoVec,speceff,"linear");
dbicmdb4 = interp1(2*sum(dbicm4,1),EsNoVec,speceff,"linear");
figure
plot(speceff,bicmdb-shannondb,altstyles{1})
hold on;
plot(speceff,ccCapdb-shannondb,altstyles{2})
plot(speceff,dbicmdb1-shannondb,altstyles{3})
plot(speceff,dbicmdb2-shannondb,altstyles{4})
plot(speceff,dbicmdb3-shannondb,altstyles{5})
plot(speceff,dbicmdb4-shannondb,altstyles{6})
grid
xlim([0.5 5])
legend('BICM','CM capacity','DBICM [0,0,1,0,0,1]','DBICM [0,1,0,0,1,0]','DBICM [1,0,1,1,0,1]','DBICM [0,1,1,0,1,1]')
ylabel('Gap to Channel Capacity (dB)')
xlabel('Spectral Efficiency (bits/s/Hz)')

function out_data = add_awgn(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end

function sumprob = prob(b,ii,y,ebno,symbollist,symlabels)

    noise_var = 1/(10^((ebno+10*log10(3))/10));
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