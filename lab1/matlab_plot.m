% Step 1: Read the data from the text file
addpath('C:\Users\User\Desktop\4061\lab1');
data = load(['lab1_bcc.txt']); % Assumes coordinates.txt is in the same directory

% Step 2: Extract x and y coordinates
x = data(:, 1); % First column for x-coordinates
y = data(:, 2); % Second column for y-coordinates
z = data(:,3);


% Step 3: Plot the coordinates
figure; % Create a new figure
scatter3(x,y,z,'filled',SizeData = 1000); % Plot with circles and lines connecting the points

xlabel('X Coordinates'); % Label for x-axis
ylabel('Y Coordinates'); % Label for y-axis
zlabel('Z');
title('Plot of Coordinates'); % Title of the plot
grid on; % Add grid for better visualization