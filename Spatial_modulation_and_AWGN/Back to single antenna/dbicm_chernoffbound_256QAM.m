clear;

M = 256;
x = (0:M-1)';
symbols = qammod(x,M,UnitAveragePower=true,PlotConstellation=false);
% symbols = real(symbols) - 1i*imag(symbols);
labels = {'00001000','00001001','00001010','00001011','00001100','00001101','00001110',...
'00001111','00000000','00000001','00000010','00000011','00000100','00000101',...
'00000110','00000111','00011000','00011001','00011010','00011011','00011100',...
'00011101','00011110','00011111','00010000','00010001','00010010','00010011',...
'00010100','00010101','00010110','00010111','00101000','00101001','00101010',...
'00101011','00101100','00101101','00101110','00101111','00100000','00100001',...
'00100010','00100011','00100100','00100101','00100110','00100111','00111000',...
'00111001','00111010','00111011','00111100','00111101','00111110','00111111',...
'00110000','00110001','00110010','00110011','00110100','00110101','00110110',...
'00110111','01001000','01001001','01001010','01001011','01001100','01001101',...
'01001110','01001111','01000000','01000001','01000010','01000011','01000100',...
'01000101','01000110','01000111','01011000','01011001','01011010','01011011',...
'01011100','01011101','01011110','01011111','01010000','01010001','01010010',...
'01010011','01010100','01010101','01010110','01010111','01101000','01101001',...
'01101010','01101011','01101100','01101101','01101110','01101111','01100000',...
'01100001','01100010','01100011','01100100','01100101','01100110','01100111',...
'01111000','01111001','01111010','01111011','01111100','01111101','01111110',...
'01111111','01110000','01110001','01110010','01110011','01110100','01110101',...
'01110110','01110111','10001000','10001001','10001010','10001011','10001100',...
'10001101','10001110','10001111','10000000','10000001','10000010','10000011',...
'10000100','10000101','10000110','10000111','10011000','10011001','10011010',...
'10011011','10011100','10011101','10011110','10011111','10010000','10010001',...
'10010010','10010011','10010100','10010101','10010110','10010111','10101000',...
'10101001','10101010','10101011','10101100','10101101','10101110','10101111',...
'10100000','10100001','10100010','10100011','10100100','10100101','10100110',...
'10100111','10111000','10111001','10111010','10111011','10111100','10111101',...
'10111110','10111111','10110000','10110001','10110010','10110011','10110100',...
'10110101','10110110','10110111','11001000','11001001','11001010','11001011',...
'11001100','11001101','11001110','11001111','11000000','11000001','11000010',...
'11000011','11000100','11000101','11000110','11000111','11011000','11011001',...
'11011010','11011011','11011100','11011101','11011110','11011111','11010000',...
'11010001','11010010','11010011','11010100','11010101','11010110','11010111',...
'11101000','11101001','11101010','11101011','11101100','11101101','11101110',...
'11101111','11100000','11100001','11100010','11100011','11100100','11100101',...
'11100110','11100111','11111000','11111001','11111010','11111011','11111100',...
'11111101','11111110','11111111','11110000','11110001','11110010','11110011',...
'11110100','11110101','11110110','11110111'};

m = log2(M);
l = 0; %No apriori

dataTable = table();

EbNodb = 0; %Eb/No value in dB
% esno = 10^((EbNodb +10*log10(m))/10);
esno = 10^((5.98)/10);

for q = 0:(2^m-1)
    q

    D  = 0;

    binaryString = dec2bin(q, m); % 4 specifies the minimum number of bits
    delayscheme = str2num(binaryString')'; % Convert binary string to binary vector
    % delayscheme = [1 1 1 0];
    
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
