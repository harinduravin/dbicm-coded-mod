lengthVector = 4;  % Set the desired length
allPossibleVectors = generateVectors(lengthVector);

% Display the generated vectors
disp('All possible vectors:');
disp(allPossibleVectors);



function allVectors = generateVectors(lengthVector)
    % Initialize an empty cell array to store the generated vectors
    allVectors = cell(0);

    % Call the recursive function to generate vectors
    generate([], lengthVector);

    % Recursive function to generate vectors
    function generate(currentVector, remainingLength)
        % Base case: if the remaining length is 0, add the current vector
        if remainingLength == 0
            allVectors{end+1} = currentVector;
            return;
        end

        % Recursive calls with different values (0, 1, 2)
        generate([currentVector, 0], remainingLength - 1);
        generate([currentVector, 1], remainingLength - 1);
        generate([currentVector, 2], remainingLength - 1);
    end
end
