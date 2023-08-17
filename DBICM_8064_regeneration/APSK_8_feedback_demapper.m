function combined_array = APSK_32_feedback_demapper(complex_numbers, ebno, extrinsic, delay_profile)


    prob_0 = 1./(1+exp(-extrinsic));
    prob_1 = 1 + (-1)./(1+exp(-extrinsic));

    symbols = { 
        6.03022689e-01+0.00000000e+00j,  3.69244903e-17+6.03022689e-01j,...
       -6.03022689e-01+7.38489806e-17j, -1.10773471e-16-6.03022689e-01j,...
        9.04534034e-01+9.04534034e-01j, -9.04534034e-01+9.04534034e-01j,...
       -9.04534034e-01-9.04534034e-01j,  9.04534034e-01-9.04534034e-01j
      };

    labels = {'011', '001', '000', '010', '111', '101', '100', '110'};
    m = 3;


    llr_matrix = get_all_llr(complex_numbers, symbols, labels, ebno, delay_profile, prob_0,m);
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

function llr = get_llr_weighted(delay_profile, prob_0, complex_numbers, symbols, labels, noise_var, position)
    llr_1 = get_indices(position, '1', char(labels));
    llr_0 = get_indices(position, '0', char(labels));
    
    symbols_llr_1 = cell2mat(symbols(llr_1));
    symbols_llr_0 = cell2mat(symbols(llr_0));

    labels_llr_1 = labels(llr_1);
    labels_llr_0 = labels(llr_0);

    
    delay_positions = find(delay_profile(:) == 1);

    all_positions_1 = zeros(length(delay_positions),length(labels_llr_1)/2,2);
    all_positions_0 = zeros(length(delay_positions),length(labels_llr_0)/2,2);

    for i = 1: length(delay_positions)
        all_positions_1(i,:,1) = get_indices(delay_positions(i), '0', char(labels_llr_1));
        all_positions_0(i,:,1) = get_indices(delay_positions(i), '0', char(labels_llr_0));
        all_positions_1(i,:,2) = get_indices(delay_positions(i), '1', char(labels_llr_1));
        all_positions_0(i,:,2) = get_indices(delay_positions(i), '1', char(labels_llr_0));
    end

    init_weight_1 = ones(length(symbols_llr_1),length(complex_numbers));
    init_weight_0 = ones(length(symbols_llr_1),length(complex_numbers));

    for i = 1: length(delay_positions)
        init_weight_1(all_positions_1(i,:,1),:) = init_weight_1(all_positions_1(i,:,1),:).*repmat(prob_0(delay_positions(i),:),length(symbols_llr_1)/2,1);
        init_weight_1(all_positions_1(i,:,2),:) = init_weight_1(all_positions_1(i,:,2),:).*repmat(1-prob_0(delay_positions(i),:),length(symbols_llr_1)/2,1);

        init_weight_0(all_positions_0(i,:,1),:) = init_weight_0(all_positions_0(i,:,1),:).*repmat(prob_0(delay_positions(i),:),length(symbols_llr_0)/2,1);
        init_weight_0(all_positions_0(i,:,2),:) = init_weight_0(all_positions_0(i,:,2),:).*repmat(1-prob_0(delay_positions(i),:),length(symbols_llr_0)/2,1);
    end
    
    eucl_llr_1_real = pdist2(real(symbols_llr_1.'), real(complex_numbers));
    eucl_llr_1_imag = pdist2(imag(symbols_llr_1.'), imag(complex_numbers));
    eucl_llr_1_squared = eucl_llr_1_real.^2 + eucl_llr_1_imag.^2;
    eucl_llr_1_eucl = -eucl_llr_1_squared/(2*noise_var);

    eucl_llr_0_real = pdist2(real(symbols_llr_0.'), real(complex_numbers));
    eucl_llr_0_imag = pdist2(imag(symbols_llr_0.'), imag(complex_numbers));
    eucl_llr_0_squared = eucl_llr_0_real.^2 + eucl_llr_0_imag.^2;
    eucl_llr_0_eucl = -eucl_llr_0_squared/(2*noise_var);

    % Weighting the numbers before summing.
    % Optimized accurate LLR
    weighted_exp_1 = exp(eucl_llr_1_eucl).*init_weight_1;
    weighted_exp_0 = exp(eucl_llr_0_eucl).*init_weight_0;  
    weighted_exp_1_sum = sum(weighted_exp_1,1);
    weighted_exp_0_sum = sum(weighted_exp_0,1); 

    llr = log(weighted_exp_1_sum) - log(weighted_exp_0_sum);


end

function llrs = get_all_llr(complex_numbers, symbols, labels, ebno, delay_profile, prob_0,m)
    llrs = zeros(m, length(complex_numbers));

    noise_var = 1/((3/4)*2*10^((ebno+10*log10(m))/10));
    for i = 1:m
        if delay_profile(i) == 1
            llrs(i,:) = get_llr(complex_numbers, symbols, labels, noise_var, i);
            % llrs(i,:) = extrinsic(i,:);
        else
            % llrs(i,:) = get_llr(complex_numbers, symbols, labels, noise_var, i);
            llrs(i,:) = get_llr_weighted(delay_profile, prob_0, complex_numbers, symbols, labels, noise_var, i);
        end
    end
end


