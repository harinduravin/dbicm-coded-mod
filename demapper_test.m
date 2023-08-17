distinctRGBColors = [
    0.298, 0.447, 0.690;   % Blue
    0.333, 0.658, 0.407;   % Green
    0.768, 0.305, 0.321;   % Red
    0.505, 0.447, 0.698;   % Purple
    0.800, 0.725, 0.454;   % Brown
    0.392, 0.709, 0.803;   % Cyan
    0.993, 0.750, 0.328;   % Orange
    0.855, 0.004, 0.278;   % Maroon
    0.552, 0.827, 0.780;   % Mint
    0.992, 0.682, 0.380;   % Salmon
    0.890, 0.101, 0.109;   % Crimson
    0.169, 0.510, 0.729;   % Sky Blue
    0.584, 0.403, 0.741;   % Lavender
    0.769, 0.580, 0.263;   % Gold
    0.729, 0.729, 0.729;   % Gray
    0.529, 0.808, 0.922;   % Light Blue
    1.000, 1.000, 0.600;   % Yellow
    0.701, 0.701, 0.701;   % Silver
    0.992, 0.753, 0.753;   % Pink
    0.498, 0.498, 0.498;   % Dark Gray
    0.737, 0.502, 0.741;   % Violet
    0.769, 0.376, 0.180;   % Rust
    0.937, 0.937, 0.937;   % Light Gray
    0.275, 0.941, 0.941;   % Aqua
    0.741, 0.467, 0.702;   % Orchid
    0.365, 0.635, 0.655;   % Teal
    0.941, 0.503, 0.749;   % Rose
    0.278, 0.439, 0.753;   % Royal Blue
    0.780, 0.780, 0.780;   % Gray
    0.620, 0.282, 0.278;   % Brick Red
    0.110, 0.337, 0.663;   % Navy Blue
    0.635, 0.078, 0.184;   % Dark Red
];

% Define the number of complex numbers you want to generate
numComplexNumbers = 10000;

% Generate random angles between 0 and 2*pi
angles = 2 * pi * rand(1, numComplexNumbers);

% Generate random radii between 0 and 1
radii = 1.9*sqrt(rand(1, numComplexNumbers));

% Create complex numbers using polar coordinates
complexNumbers = radii .* exp(1i * angles);
complexNumbers = complexNumbers.';

extrinsic = zeros(5,length(complexNumbers));

extrinsic(3,:) = 0;
extrinsic(5,:) = Inf;

% demapped = APSK_32_demapper_optimized(complexNumbers, 5);
demapped = APSK_32_feedback_demapper(complexNumbers, 5, extrinsic, [0 0 0 0 1]);
% demapped = APSK_8_demapper_optimized(complexNumbers, 5);

Harddemapped = demapped < 0;

% Convert binary values to cell array of strings
binaryString = num2str(Harddemapped);
binaryString = reshape(binaryString, 5, []).'; % Reshape into groups of 5

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


