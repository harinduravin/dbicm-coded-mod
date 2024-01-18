%16-QAM

symbols = [-0.948683298050514 +   0.948683298050514i,...
         -0.948683298050514 +   0.316227766016838i,...
         -0.948683298050514 -   0.316227766016838i,...
         -0.948683298050514 -   0.948683298050514i,...
         -0.316227766016838 +   0.948683298050514i,...
         -0.316227766016838 +   0.316227766016838i,...
         -0.316227766016838 -   0.316227766016838i,...
         -0.316227766016838 -   0.948683298050514i,...
          0.316227766016838 +   0.948683298050514i,...
          0.316227766016838 +   0.316227766016838i,...
          0.316227766016838 -   0.316227766016838i,...
          0.316227766016838 -   0.948683298050514i,...
          0.948683298050514 +   0.948683298050514i,...
          0.948683298050514 +   0.316227766016838i,...
          0.948683298050514 -   0.316227766016838i,...
          0.948683298050514 -   0.948683298050514i,...
          ];

% Define the labels
% Gray
labels = {'0000','0001','0011','0010','0100','0101','0111','0110','1100','1101','1111','1110','1000','1001','1011','1010'};

% awgn optimized
% labels = {'1101','0111','1010','0000','1011','0001','1100','0110','0100','0010','1111','1001','1110','1000','0101','0011'};

m = 4;
l = m-1; %Ideal apriori
D  = 0;
EbNodb = 0; %Eb/No value in dB
% esno = 10^((EbNodb +10*log10(m))/10);
esno = 10^((8)/10);

for i = 1:m
    for b = 0:1
        symbol_subset = symbols(get_symbol_indices(symbols, labels, i, b));
        for t = 1:length(symbol_subset)

            symbol_subset_bar = symbols(get_symbol_indices(symbols, labels, i, ~b));
            label_subset_bar = labels(get_symbol_indices(symbols, labels, i, ~b)); 
            linearIndex = find(symbols==symbol_subset(t));
            labelitem = char(labels(linearIndex));

            aprioripositions = generateVector(i, m, l);
            aprioribits = labelitem(logical(aprioripositions));
            indexlist = get_dbicmindices(aprioripositions, aprioribits, char(label_subset_bar));
            symbol_subset_bar = symbol_subset_bar(indexlist);

            for s = 1:length(symbol_subset_bar)
                length(symbol_subset_bar)
                D = D + (1/(m*2^m))*(exp(-(esno/4)*abs(symbol_subset(t)-symbol_subset_bar(s))^2));
            end
        end
    end
end

D

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
