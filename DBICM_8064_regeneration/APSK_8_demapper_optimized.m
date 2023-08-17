function combined_array = APSK_8_demapper_optimized(complex_numbers, ebno)
    
    symbols = { 
        6.03022689e-01+0.00000000e+00j,  3.69244903e-17+6.03022689e-01j,...
       -6.03022689e-01+7.38489806e-17j, -1.10773471e-16-6.03022689e-01j,...
        9.04534034e-01+9.04534034e-01j, -9.04534034e-01+9.04534034e-01j,...
       -9.04534034e-01-9.04534034e-01j,  9.04534034e-01-9.04534034e-01j
      };

    labels = {'011', '001', '000', '010', '111', '101', '100', '110'};
    m = 3;

    llr_matrix = get_all_llr(complex_numbers, symbols, labels, ebno, m);
    combined_array =  -llr_matrix(:);


end

function indices = get_indices(position, bit_value, labels)
    indices = find(labels(:, position) == bit_value);
end

function llr = get_llr(complex_numbers, symbols, labels, noise_var, position)
    llr_1 = get_indices(position, '1', char(labels));
    llr_0 = get_indices(position, '0', char(labels));
    
    symbols_llr_1 = cell2mat(symbols(llr_1));
    symbols_llr_0 = cell2mat(symbols(llr_0));
    
    eucl_llr_1_real = pdist2(real(symbols_llr_1.'), real(complex_numbers));
    eucl_llr_1_imag = pdist2(imag(symbols_llr_1.'), imag(complex_numbers));
    eucl_llr_1_squared = eucl_llr_1_real.^2 + eucl_llr_1_imag.^2;
    eucl_llr_1_eucl = -eucl_llr_1_squared/(2*noise_var);

    eucl_llr_0_real = pdist2(real(symbols_llr_0.'), real(complex_numbers));
    eucl_llr_0_imag = pdist2(imag(symbols_llr_0.'), imag(complex_numbers));
    eucl_llr_0_squared = eucl_llr_0_real.^2 + eucl_llr_0_imag.^2;
    eucl_llr_0_eucl = -eucl_llr_0_squared/(2*noise_var);

    % Optimized accurate LLR
    llr = log(sum(exp(eucl_llr_1_eucl),1)) - log(sum(exp(eucl_llr_0_eucl),1));
end

function llrs = get_all_llr(complex_numbers, symbols, labels, ebno, m)
    llrs = zeros(m, length(complex_numbers));

    noise_var = 1/((3/4)*2*10^((ebno+10*log10(m))/10));
    for i = 1:m
        llrs(i,:) = get_llr(complex_numbers, symbols, labels, noise_var, i);
    end
end


