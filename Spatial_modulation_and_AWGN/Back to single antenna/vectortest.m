iValue = 5;
mValue = 5;
lValue = 0;
resultVector = generateVector(iValue, mValue, lValue);
disp(resultVector);


function resultVector = generateVector(i, m, l)

    % if ~(1 <= i && i <= m && 1 <= l && l <= m)
    %     error('Invalid values for i, m, or l');
    % end

    resultVector = zeros(1, m);  % Initialize vector with all 1s
    resultVector(i) = 0;        % Set the ith element to 0
    filtervec = resultVector([1:i-1,i+1:end]);
    filtervec(1:l) = 1;      % Set elements from 1st to lth to 1
    resultVector([1:i-1,i+1:end]) = filtervec;
end


