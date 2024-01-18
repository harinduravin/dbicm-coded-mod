clear;
EbNoVec = (-15:2:20)';      % Eb/No values (dB)
testval = 0;
testval2 = 0;
N = 90000;



% 2x2 mimo channel
Nr = 4;
Nt = 4;





% H = 1/sqrt(2)*(ones(Nr,Nt) + 1i*ones(Nr,Nt));


% Number of selected antennas from Tx antennas
Na = 2;
% bits per symbol
m = 3;
p_1 = m; % 8 PSK
% Spatial bits
p_11 = floor(log2(nchoosek(Nt,Na)));
% Total bits per hypersymbol
hsymbolbits = p_11 + Na*p_1;

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
hsymbmsg = zeros(Nt,N);
for q = 1:N
    index = find(strcmp(hsymlabels, sprintf('%d',msg(:,q)')));
    hsymbmsg(:,q) = hsymbollist(:,index);
end

capEst = zeros(hsymbolbits,length(EbNoVec));

for ii = 1:hsymbolbits
    for n = 1:length(EbNoVec)
        EbNoVec(n)
    
        llccap = 0;
       
        % Reset the error and bit counters
        numErrs = 0;
        numBits = 0;
    
        % accumulate=zeros(2,2,N/2);

        tic

        parfor j = 1:(N)
    
            H = 1/sqrt(2)*(randn(Nr,Nt) + 1i*randn(Nr,Nt));
    
            symb = hsymbmsg(:,j);
            
            r = H*hsymbmsg(:,j);
    
            % AWGN noise at receiver
            y = add_awgn(r, EbNoVec(n),1,1);
    
            % Bit levels 
            bitlevels = msg(:,j);
            if (ii == 3)&&(n == 9)
                testval = testval + (1/N)*prob(bitlevels(ii),ii,y,EbNoVec(n),m,H,hsymbollist,hsymlabels);
            end
            if (ii == 6)&&(n == 9)
                testval2 = testval2 + (1/N)*prob(bitlevels(ii),ii,y,EbNoVec(n),m,H,hsymbollist,hsymlabels);
            end    
            llccap = llccap + (1/N)*(1-log2(prob(bitlevels(ii),0,y,EbNoVec(n),m,H,hsymbollist,hsymlabels)/prob(bitlevels(ii),ii,y,EbNoVec(n),m,H,hsymbollist,hsymlabels)));
        end
        toc
        % Store the capacity
        capEst(ii,n) = llccap;
    end
end
altstyles = {'bo--'; 'rs-.'; 'kd-'; 'g+--'; 'kv-'; 'b--'; 'r-'; 'k-'; 'g--'; 'g-';'b-'; 'mv-'};
Legend=cell(hsymbolbits,1);
for iter=1:hsymbolbits
    Legend{iter}=strcat('Bit -', num2str(iter));
end

for z = 1:hsymbolbits
    plot(EbNoVec,capEst(z,:),altstyles{z})
    hold on;
end
grid
legend(Legend)
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

function sumprob = prob(b,ii,y,ebno,m,H,hsymbollist,hsymlabels)

    noise_var = 1/(2*10^((ebno+10*log10(1))/10));
    sumprob = get_sum(y, hsymbollist, hsymlabels, noise_var, ii,b,H,m);

end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
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

