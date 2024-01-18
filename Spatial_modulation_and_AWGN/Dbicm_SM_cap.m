clear;
EbNoVec = (-15:2:20)';      % Eb/No values (dB)
testval = 0;
testval2 = 0;

N_bicm = 32000;
N = N_bicm*16/256;


% 4x4 mimo channel
Nr = 4;
Nt = 4;

% Number of selected antennas from Tx antennas
Na = 2;
% bits per symbol
m = 3;
p_1 = m; % 8 PSK
% Spatial bits
p_11 = floor(log2(nchoosek(Nt,Na)));
% Total bits per hypersymbol
hsymbolbits = p_11 + Na*p_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delayscheme = zeros(1,hsymbolbits);

delayscheme(3) = 1;
delayscheme(4) = 1;
delayscheme(6) = 1;
delayscheme(7) = 1;

numDelaybits = sum(delayscheme(:) == 1);

delaylabels = cell(1, 2^numDelaybits);

for i = 0:(2^numDelaybits - 1)
    binary_seq = dec2bin(i, numDelaybits);
    delaylabels{i + 1} = binary_seq;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hsymlabels = cell(1, 2^hsymbolbits);

for i = 0:(2^hsymbolbits - 1)
    binary_seq = dec2bin(i, hsymbolbits);
    hsymlabels{i + 1} = binary_seq;
end

binary_matrix = zeros(hsymbolbits, 2^hsymbolbits);
for i = 1:2^hsymbolbits
    hsymlabel = hsymlabels{i};
    for j = 1:hsymbolbits
        binary_matrix(j, i) = str2double(hsymlabel(j));
    end
end

hsymbollist = zeros(Nt,2^hsymbolbits);
for j = 1:2^hsymbolbits
    % Bit levels considered
    allbits = int8(binary_matrix(:,j)');
    
    spatialbits = allbits(1:p_11);
    bitlevels = reshape(allbits(p_11+1:end), [], 1);
    
    % Transmit symbols for 2x2 mimo
    x_na = pskmod(bitlevels,8,0,'gray','InputType','bit');
    
    % Spatial modulation for two bit spatial modulation
    if all(spatialbits == [0 0])
        hsymbollist(:,j) = [0 1;1 0;0 0;0 0]*x_na;
    elseif all(spatialbits == [0 1])
        hsymbollist(:,j) = [0 1;0 0;1 0;0 0]*x_na;
    elseif all(spatialbits == [1 0])
        hsymbollist(:,j) = [0 0;0 1;0 0;1 0]*x_na;
    elseif all(spatialbits == [1 1])
        hsymbollist(:,j) = [0 0;0 0;0 1;1 0]*x_na;
    end
end

msg = logical(randi([0 1],hsymbolbits,N));

hsymbmsg = zeros(Nt,N,2^numDelaybits);
for i = 1:2^numDelaybits
    temporary_binary_matrix = msg;
    temporary_binary_matrix(logical(delayscheme),:) = repmat(str2num(char(delaylabels(i))'),1,N);
    for q = 1:N
        index = find(strcmp(hsymlabels, sprintf('%d',temporary_binary_matrix(:,q)')));
        hsymbmsg(:,q,i) = hsymbollist(:,index);
    end
end

msg_bicm = logical(randi([0 1],hsymbolbits,N_bicm));
hsymbmsgbicm = zeros(Nt,N_bicm);
for q = 1:N_bicm
    index = find(strcmp(hsymlabels, sprintf('%d',msg_bicm(:,q)')));
    hsymbmsgbicm(:,q) = hsymbollist(:,index);
end

capEst = zeros(hsymbolbits,length(EbNoVec));

for ii = 1:hsymbolbits
    for n = 1:length(EbNoVec)
        EbNoVec(n)
    
        llccap = 0;

        if delayscheme(ii) == 1
            for j = 1:N_bicm
        
                H = 1/sqrt(2)*(randn(Nr,Nt) + 1i*randn(Nr,Nt));
        
                % symb = hsymbmsg(:,j);
                
                r = H*hsymbmsgbicm(:,j);
        
                % AWGN noise at receiver
                y = add_awgn(r, EbNoVec(n),1,1/2);
        
                % Bit levels 
                bitlevels = msg_bicm(:,j);

                llccap = llccap + (1/N_bicm)*(log2(prob(bitlevels(ii),0,y,EbNoVec(n),m,H,hsymbollist,hsymlabels)/prob(bitlevels(ii),ii,y,EbNoVec(n),m,H,hsymbollist,hsymlabels)));
            end
        else
            for zz = 1:length(delaylabels)
                delayedbits = delaylabels(zz);
                indexlist = get_dbicmindices(delayscheme, char(delayedbits), char(hsymlabels));
        
                delayhsymlabels = hsymlabels(indexlist);
                delayhsymbollist = hsymbollist(:,indexlist);
    
                cap = 0;
    
                for j = 1:(N)
            
                    H = 1/sqrt(2)*(randn(Nr,Nt) + 1i*randn(Nr,Nt));
            
                    % symb = hsymbmsg(:,j);
        
                    r = H*hsymbmsg(:,j,zz);
        
                    % AWGN noise at receiver
                    y = add_awgn(r, EbNoVec(n),1,1/2);
        
                    % Bit levels
                    bitlevels = msg(:,j);
        
                    cap = cap + (1/N)*(log2(prob(bitlevels(ii),0,y,EbNoVec(n),m,H,delayhsymbollist,delayhsymlabels)/prob(bitlevels(ii),ii,y,EbNoVec(n),m,H,delayhsymbollist,delayhsymlabels)));
                end
                llccap = llccap + (1/length(delaylabels))*cap;
            end
            
        end
    
        % Store the capacity
        capEst(ii,n) = 1-llccap;
    end
end
altstyles = {'bo--'; 'rs-.'; 'kd-'; 'g+--'; 'kv-'; 'b--'; 'r-'; 'k-'; 'g--'; 'g-';'b-'; 'mv-'};
Legend=cell(hsymbolbits,1);
for iter=1:hsymbolbits
    Legend{iter}=strcat('Bit -', num2str(iter));
end

figure
for z = 1:hsymbolbits
    plot(EbNoVec,capEst(z,:),altstyles{z})
    hold on;
end
grid
legend(Legend)
xlabel('Es/No (dB)')
ylabel('Bit capacity')

figure
plot(EbNoVec,sum(capEst,1),altstyles{9})
grid
legend('[0,0,1,0,1,1,0,1]')
xlabel('Es/No (dB)')
ylabel('Total capacity')

function out_data = add_awgn(signal, ebno_db, code_rate,bits)
    % Calculate the noise power
    ebno = 10^((ebno_db+10*log10(bits))/10);
    noise_pow = 1/sqrt(2*ebno);
    
    % Generate complex Gaussian noise
    noise = noise_pow * (1/sqrt(code_rate)) * (randn(size(signal)) + 1i * randn(size(signal)));
    
    % Add noise to the signal
    out_data = signal + noise;
end

function sumprob = prob(b,ii,y,ebno,m,H,hsymbollist,hsymlabels)

    noise_var = 1/(2*10^((ebno+10*log10(1/2))/10));
    sumprob = get_sum(y, hsymbollist, hsymlabels, noise_var, ii,b,H,m);

end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
end

function indices = get_dbicmindices(delayscheme, delayedbits, labels)
    indices = find(all(labels(:, logical(delayscheme)) == delayedbits,2));
end

function sumval = get_sum(r, symbols, labels, noise_var, position, value,H,m)

    p_11 = 1; %number of spatial bits
    
    
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
    % sel_hyp_symbols = [Array1(:), Array2(:)];
    Ha = H*selected_symbols;
    r_Ha = r -Ha;
    eucl_norms = vecnorm(r_Ha).^2;
    sumval = sum(exp(-eucl_norms/(2*noise_var)));

    % eucl_norms = eucl_dist_norm(Ha, r, noise_var);
    % 
    % sumval = sum(exp(eucl_norms));
end

