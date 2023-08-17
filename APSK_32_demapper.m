function combined_array = APSK_32_demapper(complex_numbers, ebno)
    
    symbols = { 
        0.164079178283360 + 0.164079178283360i,...
        -0.164079178283360 + 0.164079178283360i,...
        -0.164079178283360 - 0.164079178283360i,...
        0.164079178283360 - 0.164079178283360i,...
        0.638788528436365 + 0.171162870328789i,...
        0.467625658107576 + 0.467625658107576i,...
        0.171162870328789 + 0.638788528436365i,...
        -0.171162870328789 + 0.638788528436365i,...
        -0.467625658107576 + 0.467625658107576i,...
        -0.638788528436365 + 0.171162870328789i,...
        -0.638788528436365 - 0.171162870328789i,...
        -0.467625658107576 - 0.467625658107576i,...
        -0.171162870328790 - 0.638788528436365i,...
        0.171162870328789 - 0.638788528436365i,...
        0.467625658107576 - 0.467625658107576i,...
        0.638788528436365 - 0.171162870328790i,...
        1.26309318727039 + 0.251244856101071i,...
        1.07079869947673 + 0.715484816502116i,...
        0.715484816502116 + 1.07079869947673i,...
        0.251244856101071 + 1.26309318727039i,...
        -0.251244856101070 + 1.26309318727039i,...
        -0.715484816502115 + 1.07079869947673i,...
        -1.07079869947673 + 0.715484816502116i,...
        -1.26309318727039 + 0.251244856101071i,...
        -1.26309318727039 - 0.251244856101071i,...
        -1.07079869947673 - 0.715484816502115i,...
        -0.715484816502116 - 1.07079869947673i,...
        -0.251244856101071 - 1.26309318727039i,...
        0.251244856101071 - 1.26309318727039i,...
        0.715484816502116 - 1.07079869947673i,...
        1.07079869947673 - 0.715484816502116i,...
        1.26309318727039 - 0.251244856101071i
      };

    labels = {'01111', '01101', '11101', '11111', '01110', '00110', '00111',...
               '00101', '00100', '01100', '11100', '10100', '10101', '10111',...
               '10110', '11110', '01011', '01010', '00010', '00011', '00001',...
               '00000', '01000', '01001', '11001', '11000', '10000', '10001',...
               '10011', '10010', '11010', '11011'};

    combined_array = zeros(length(complex_numbers)*5,1);
    for i = 1:length(complex_numbers)
        abc = get_all_llr(complex_numbers(i), symbols, labels, ebno);
        combined_array(5*i-4:5*i) = -abc;
    end

end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
end

function eucl_dist = eucl_dist_norm(z1, z2, noise_var)
    real_diff = real(z1) - real(z2);
    imag_diff = imag(z1) - imag(z2);
    squared_diff = real_diff.^2 + imag_diff.^2;
    eucl_dist = -squared_diff / noise_var;
end

function llr = get_llr(complex_number, symbols, labels, noise_var, position)
    llr_1 = get_indices(position, '1', char(labels));
    llr_0 = get_indices(position, '0', char(labels));
    
    symbols_llr_1 = cell2mat(symbols(llr_1));
    symbols_llr_0 = cell2mat(symbols(llr_0));
    
    eucl_llr_1 = eucl_dist_norm(symbols_llr_1, complex_number, noise_var);
    eucl_llr_0 = eucl_dist_norm(symbols_llr_0, complex_number, noise_var);
    
    % Max log approximation
    % llr = max(eucl_llr_1) - max(eucl_llr_0);
    
    % Accurate calculation
    llr = log(sum(exp(eucl_llr_1))) - log(sum(exp(eucl_llr_0)));
end

function llrs = get_all_llr(complex_number, symbols, labels, ebno)
    llrs = zeros(5, 1);

    noise_var = 1/(2*10^((ebno+10*log10(5))/10));
    for i = 1:5
        llrs(i) = get_llr(complex_number, symbols, labels, noise_var, i);
    end
end


