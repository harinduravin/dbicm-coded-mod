clear
M = 16;
x = (0:M-1)';
symbols = dvbsapskmod(x,M,'s2x','90/180',UnitAveragePower=true,PlotConstellation=false);

% Define the labels
% Gray labelled APSK-16
labels = {'0000','1000','0100','1100','0010','1010','0110','1110','0001','1001','0101','1101','0011','1011','0111','1111'};
% 
% Plot the complex numbers
figure;
hold on;
plot(real(symbols), imag(symbols), 'o'); % Plot real vs imaginary part
text(real(symbols), imag(symbols), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','right'); % Add labels

% Set axis labels
xlabel('Real');
ylabel('Imaginary');
title('Complex Numbers with Labels');
ylim([-1.5 1.5])
xlim([-1.5 1.5])
% Set aspect ratio to be equal
axis equal;

% Display grid
grid on;

% Show plot
hold off;

% demod = dvbsapskdemod(symbols,M,'s2x','90/180',UnitAveragePower=true,OutputType='bit');
% bits = reshape(demod, [log2(M), M]);
% cellArrayOfStrings = cellstr(num2str(bits'));
