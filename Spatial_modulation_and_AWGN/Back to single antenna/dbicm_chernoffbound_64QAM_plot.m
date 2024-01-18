% clear;

load('64qam_gap.mat','dbicmdb1','dbicmdb2','dbicmdb3','dbicmdb4','speceff');
% 'DBICM [0,0,1,0,0,1]','DBICM [0,1,0,0,1,0]','DBICM [1,0,1,1,0,1]','DBICM [0,1,1,0,1,1]'

M = 64;
x = (0:M-1)';
symbols = qammod(x,M,UnitAveragePower=true,PlotConstellation=false);
symbols = real(symbols) - 1i*imag(symbols);
labels = {'000000','000001','000010','000011','000100','000101','000110','000111','001000','001001','001010','001011','001100',...
    '001101','001110','001111','010000','010001','010010','010011','010100','010101','010110','010111','011000','011001','011010',...
    '011011','011100','011101','011110','011111','100000','100001','100010','100011','100100','100101','100110','100111','101000',...
    '101001','101010','101011','101100','101101','101110','101111','110000','110001','110010','110011','110100','110101','110110',...
    '110111','111000','111001','111010','111011','111100','111101','111110','111111'};

m = 6;
l = 0; %No apriori

% dataTable = table();

% EbNodb = 0; %Eb/No value in dB
% esno = 10^((EbNodb +10*log10(m))/10);
% esno = 10^((5)/10);
metrics = zeros(4,length(dbicmdb1));
delayschemes = [0,0,1,0,0,1;0,1,0,0,1,0;1,0,1,1,0,1;0,1,1,0,1,1];
dbicmdbs = zeros(4,length(dbicmdb1));
dbicmdbs(1,:) = dbicmdb1;
dbicmdbs(2,:) = dbicmdb2;
dbicmdbs(3,:) = dbicmdb3;
dbicmdbs(4,:) = dbicmdb4;

for z = 1:4

    delayscheme = delayschemes(z,:);
    
    for q = 1:length(dbicmdb1)
    
        D  = 0;
        esno = 10^((dbicmdbs(z,q))/10);
        
        % binaryString = dec2bin(q, m); % 4 specifies the minimum number of bits
        % delayscheme = str2num(binaryString')'; % Convert binary string to binary vector
    
        
        bitwisecosts = zeros(size(delayscheme));
        
        for i = 1:m
            Dm = 0;
            for b = 0:1
                symbol_subset = symbols(get_symbol_indices(symbols, labels, i, b));
                for t = 1:length(symbol_subset)
        
                    symbol_subset_bar = symbols(get_symbol_indices(symbols, labels, i, ~b));
                    label_subset_bar = labels(get_symbol_indices(symbols, labels, i, ~b)); 
                    linearIndex = find(symbols==symbol_subset(t));
                    labelitem = char(labels(linearIndex));
        
                    if delayscheme(i)>0
                        aprioripositions = generateVector(i, m, l);
                        aprioribits = labelitem(logical(aprioripositions));
                        indexlist = get_dbicmindices(aprioripositions, aprioribits, char(label_subset_bar));
                        symbol_subset_bar = symbol_subset_bar(indexlist);
                        % label_subset_bar = label_subset_bar(indexlist);
                    else
                        aprioripositions = delayscheme;
                        aprioribits = labelitem(logical(aprioripositions));
                        indexlist = get_dbicmindices(aprioripositions, aprioribits, char(label_subset_bar));
                        symbol_subset_bar = symbol_subset_bar(indexlist);
                    end
    
                    % Dm = Dm + (1/(m*2^m))*qfunc(abs(symbol_subset(t)-nearest_neighbour(symbol_subset(t),symbol_subset_bar))/sqrt(2*esno));
                    Dm = Dm + (1/(m*2^m))*(exp(-(esno/4)*abs(symbol_subset(t)-nearest_neighbour(symbol_subset(t),symbol_subset_bar))^2));
                    % Dm = Dm + (1/(m*2^m))*(1/abs(symbol_subset(t)-nearest_neighbour(symbol_subset(t),symbol_subset_bar))^2);
                    
                    % for s = 1:length(symbol_subset_bar)
                    %     % Dm = Dm +
                    %     % (1/(m*2^m))*(exp(-(esno/4)*abs(symbol_subset(t)-symbol_subset_bar(s))^2));
                    %     % AWGN
                    %     Dm = Dm + (1/(m*2^m))*(1/abs(symbol_subset(t)-symbol_subset_bar(s))^2); % HMMSED
                    % end
                end
            end
            bitwisecosts(i) = Dm;
            D = D + Dm;
        end
    
        metrics(z,q) = D;
        
        % delayscheme
        % bitwisecosts
        % D
        % Generate data for a new row (replace this with your data generation logic)
        % newRow = table({sprintf('%d', delayscheme)}, D, bitwisecosts(1), bitwisecosts(2), bitwisecosts(3), bitwisecosts(4) ,'VariableNames', {'Delay Scheme', 'Overall', '1', '2', '3', '4'});
        
        % Add the new row to the data table
        % dataTable = [dataTable; newRow];
    end
end

altstyles = {'bo--'; 'rs-.'; 'kd-'; 'g+--'; 'kv-'; 'b--'; 'r-'; 'k-'; 'g--'; 'g-';'b-'; 'mv-'};

figure
for i = 1:4
plot(speceff,metrics(i,:),altstyles{i+2})
hold on;
end
grid
legend('DBICM [0,0,1,0,0,1]','DBICM [0,1,0,0,1,0]','DBICM [1,0,1,1,0,1]','DBICM [0,1,1,0,1,1]')
% xlim([0.5 5])
% legend('BICM','CM capacity','DBICM [0,0,1,0,0,1]','DBICM [0,1,0,0,1,0]','DBICM [1,0,1,1,0,1]','DBICM [0,1,1,0,1,1]')
ylabel('metric cost')
xlabel('Spectral Efficiency (bits/s/Hz)')

% sortedtable = sortrows(dataTable,'Overall','ascend');

function symbol_indices = get_symbol_indices(symbols, labels, position, value)
    if position == 0
        symbol_indices = 1:length(symbols);
        % selected_symbols = symbols;
    else
        if value == 0
            valuestr = '0';
        else
            valuestr = '1';
        end
        symbol_indices = get_indices(position, valuestr, char(labels));
        % selected_symbols = symbols(symbol_indices);
    end
end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
end

function indices = get_dbicmindices(delayscheme, delayedbits, labels)
    indices = find(all(labels(:, logical(delayscheme)) == delayedbits,2));
end

function resultVector = generateVector(i, m, l)

    resultVector = zeros(1, m);  % Initialize vector with all 1s
    resultVector(i) = 0;        % Set the ith element to 0
    filtervec = resultVector([1:i-1,i+1:end]);
    filtervec(1:l) = 1;      % Set elements from 1st to lth to 1
    resultVector([1:i-1,i+1:end]) = filtervec;
end

function neighbour = nearest_neighbour(x,symbols)
    % Compute the Euclidean distances
    distances = abs(symbols - x);

    % Find the index of the complex number with the minimum distance
    [~, index] = min(distances);

    % Retrieve the complex number with the minimum distance
    neighbour = symbols(index);
end
