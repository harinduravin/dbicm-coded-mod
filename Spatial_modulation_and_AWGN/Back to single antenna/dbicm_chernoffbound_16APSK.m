clear;

M = 16;
x = (0:M-1)';
symbols = dvbsapskmod(x,M,'s2x','90/180',UnitAveragePower=true,PlotConstellation=false);

% Gray labelled APSK-16
labels = {'0000','1000','0100','1100','0010','1010','0110','1110','0001','1001','0101','1101','0011','1011','0111','1111'};

m = log2(M);
l = 0; %No apriori

dataTable = table();

% EbNodb = 0; %Eb/No value in dB
% esno = 10^((EbNodb +10*log10(m))/10);
esno = 10^((0.6)/10);

for q = 0:(2^m-1)

    D  = 0;

    binaryString = dec2bin(q, m); % 4 specifies the minimum number of bits
    delayscheme = str2num(binaryString')'; % Convert binary string to binary vector
    
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
    
    % delayscheme
    % bitwisecosts
    % D
    % Generate data for a new row (replace this with your data generation logic)
    newRow = table({sprintf('%d', delayscheme)}, D, bitwisecosts(1), bitwisecosts(2), bitwisecosts(3), bitwisecosts(4) ,'VariableNames', {'Delay Scheme', 'Overall', '1', '2', '3', '4'});
    
    % Add the new row to the data table
    dataTable = [dataTable; newRow];
end

sortedtable = sortrows(dataTable,'Overall','ascend');

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
