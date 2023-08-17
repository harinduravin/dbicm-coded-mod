distinctRGBColors = [
    0.298, 0.447, 0.690;   % Blue
    0.333, 0.658, 0.407;   % Green
    0.768, 0.305, 0.321;   % Red
    0.505, 0.447, 0.698;   % Purple
    0.800, 0.725, 0.454;   % Brown
    0.392, 0.709, 0.803;   % Cyan
    0.993, 0.750, 0.328;   % Orange
];

% Define the number of complex numbers you want to generate
numComplexNumbers = 10000;

% Generate random angles between 0 and 2*pi
angles = 2 * pi * rand(1, numComplexNumbers);

% Generate random radii between 0 and 1
radii = 0.6*sqrt(rand(1, numComplexNumbers));

% Create complex numbers using polar coordinates
complexNumbers = radii .* exp(1i * angles);
complexNumbers = complexNumbers.';

% extrinsic = zeros(5,length(complexNumbers));
% 
% extrinsic(3,:) = 0;
% extrinsic(5,:) = Inf;

% demapped = APSK_32_demapper_optimized(complexNumbers, 5);
% demapped = APSK_32_feedback_demapper(complexNumbers, 5, extrinsic, [0 0 0 0 1]);
demapped = APSK_8_demapper_optimized(complexNumbers, 5);

Harddemapped = demapped < 0;

% Convert binary values to cell array of strings
binaryString = num2str(Harddemapped);
binaryString = reshape(binaryString, 3, []).'; % Reshape into groups of 3

% Convert cell array of strings to decimal values
decimalValues = bin2dec(cellstr(binaryString))';



colour_matrix = zeros(length(decimalValues),3);

colour_count = 0;

for i = 1:length(distinctRGBColors)
    colour_matrix(decimalValues == i,:) = repmat(distinctRGBColors(i,:),length(colour_matrix(decimalValues == i,:)),1);
    if length(colour_matrix(decimalValues == i,:)) > 10
        colour_count = colour_count + 1;
        length(colour_matrix(decimalValues == i,:))
    end
end

% Plot the complex numbers
figure;
scatter(real(complexNumbers), imag(complexNumbers),5,colour_matrix, 'filled');
xlabel('Real Part');
ylabel('Imaginary Part');
title(strcat('Uniformly Distributed Complex Numbers with  ',string(colour_count)));
axis equal;
grid on;


