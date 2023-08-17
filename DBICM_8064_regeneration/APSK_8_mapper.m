function output_symbols = APSK_8_mapper(binary_sequence)

    m = 3;

    complex_numbers = { 
        6.03022689e-01+0.00000000e+00j,  3.69244903e-17+6.03022689e-01j,...
       -6.03022689e-01+7.38489806e-17j, -1.10773471e-16-6.03022689e-01j,...
        9.04534034e-01+9.04534034e-01j, -9.04534034e-01+9.04534034e-01j,...
       -9.04534034e-01-9.04534034e-01j,  9.04534034e-01-9.04534034e-01j
      };

    labels = {'011', '001', '000', '010', '111', '101', '100', '110'};

    if mod(length(binary_sequence), m) ~= 0
        error('The binary sequence length must be a multiple of m for 2^(m)-APSK conversion.');
    end
    
    % Create a map object to store bit strings as keys and their corresponding symbols as values
    bit_to_symbol_map = containers.Map(labels, complex_numbers);
    
    % Convert the input bit stream to a matrix of groups of m bits
    bit_groups = reshape(binary_sequence, m, []).';
    
    % Convert each group of bits to its corresponding symbol using the mapping
    output_symbols_flat = zeros(1, size(bit_groups, 1));
    for i = 1:size(bit_groups, 1)
        bit = num2str(bit_groups(i, :));
        
        output_symbols_flat(i) = bit_to_symbol_map(strrep(bit, ' ', ''));
    end

    output_symbols = output_symbols_flat.';
end