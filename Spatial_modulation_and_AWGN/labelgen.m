n = 4;
labels = cell(1, 2^n);

for i = 0:(2^n - 1)
    binary_seq = dec2bin(i, n);
    labels{i + 1} = binary_seq;
end

% Display the labels
for i = 1:length(labels)
    fprintf('Label: %s\n', labels{i});
end


