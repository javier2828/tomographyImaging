%Javier Salazar - Radon Transform - RadonTransform('phantom.jpg', 1)
function [] = RadonTransform(image, thetaStep) % default above will use 180 angles on image. Does not return anything, must use Phantom(...) FIRST.
tic % start timer
angleArray = 0:thetaStep:179; % all angles that it will be tested at
grid = imread(image); % transform image to matrix grid
radonImage = radon(grid, angleArray);
[xRow,yRow] = size(grid); % describes dimensions of grid columns and rows
diagPixels = sqrt(xRow^2 + yRow^2); % length of the diagnal of the grid, needed to pad image so that information isn't lost due to crop
xPadding = ceil(diagPixels - xRow + 2); % amount of padding we need to add for columns
yPadding = ceil(diagPixels - yRow + 2); % amount of padding we need to add for rows
paddedGrid = zeros(xRow+xPadding, yRow+yPadding); % new bigger grid defined
xPos = 1; % variables to transfer old grid to new grid
yPos = 1; % starting position of old grid
for j = ceil(xPadding/2):(ceil(xPadding/2)+xRow-1) % define starting and ending position of new grid to fill in from old grid
    for k = ceil(yPadding/2):(ceil(yPadding/2)+yRow-1) % define range of new grid position to input old values
        paddedGrid(j,k) = grid(xPos,yPos); % transfer old grid values to new padded grid
        yPos = yPos + 1; % defines what position will be filled
    end % end second for loop of rows
    yPos = 1; % reset for new column
    xPos = xPos + 1; % defines what position will be filled
end % end first for loop of columns
numRays = size(paddedGrid,1); % the rays for all rows in the grid
lengthLine = size(paddedGrid,2); % the x dimension for the matrix. 
raySpacing = 2/(numRays-1); % space between rays
vectorRays = -1:raySpacing:1; % location of the rays on a plot from -1 to 1
[posX, posY] = meshgrid(vectorRays); % square matrices created with dimensions of length(vectorRays) that contain ray location from -1 to 1
numAngles = length(angleArray); % number of angles given the thetaStep
radonGrid = zeros(lengthLine,numAngles); % final matrix dimensions based on thetastep ONLY. Number of rays BASED ON padded grid!
for j = 1:numAngles % will need to add up results for all angles
    xPrime = (cosd(180-angleArray(j))*posX - sind(180-angleArray(j))*posY); % orthogonal transformation for new x axis based on angle
	yPrime = (sind(180-angleArray(j))*posX + cosd(180-angleArray(j))*posY); % orthogonal transformation for new rotated y axis.
	interpolatedGrid = interp2(posX,posY,paddedGrid,xPrime,yPrime,'cubic'); % EXPLAIN HERE!!!
    radonGrid(:,j) = flipud(transpose(nansum(interpolatedGrid))); % EXPLAIN HERE!!!
end % end summation of values at different angles
radonGrid = radonGrid./max(radonGrid(:)); % normalize everything so that max is 1 for greyscale image
imwrite(radonGrid, 'radonALG2Phantom.jpg'); % matrix written to image file
figure('Name','Radon ALG-X Phantom','NumberTitle','off') % Title new figure window
%imshow('radonALG2Phantom.jpg','Colormap',parula(100)) % import the image with parula colormap instead of greyscale for visual apperance
imagesc(radonGrid);
title(['Radon ALG ',num2str(numAngles),' X ',num2str(lengthLine)]); % title of plot with dimensions of radon image
h = colorbar; % display colorbar with picture to indicate density
ylabel(h, 'Density'); % label colorbar
figure('Name','Radon MATLAB Phantom','NumberTitle','off') % Title new figure window
imagesc(radonImage);
%imwrite(radonImage, 'radonMATLABPhantom.jpg');
%imshow('radonMATLABPhantom.jpg','Colormap',parula(100))
title(['Radon MATLAB Function ',num2str(size(radonImage,2)),' X ',num2str(size(radonImage,1))]);
ylabel(h, 'Density'); % label colorbar
toc % end timer
end % Radon function ends


%{
function pointIndex4 = endPoints(j,k)
global xRow;
global yRow;
global rayArray;
global angleArray;
pointIndex1 = zeros(1,8);
if (rayArray(k)>=0)
    pointIndex1(2) = yRow/2;
    pointIndex1(1) = (rayArray(k)*sind(angleArray(j))*tand(angleArray(j)) + rayArray(k)*cosd(angleArray(j)) - pointIndex1(2)*tand(angleArray(j)));
    pointIndex1(3) = -xRow/2;
    pointIndex1(4) = (-cotd(angleArray(j))*pointIndex1(3) + rayArray(k)*sind(angleArray(j)) + cotd(angleArray(j))*rayArray(k)*cosd(angleArray(j)));
    pointIndex1(5) = xRow/2;
    pointIndex1(6) = (-cotd(angleArray(j))*pointIndex1(5) + rayArray(k)*sind(angleArray(j)) + cotd(angleArray(j))*rayArray(k)*cosd(angleArray(j)));
    pointIndex1(8) = -yRow/2;
    pointIndex1(7) = (rayArray(k)*sind(angleArray(j))*tand(angleArray(j)) + rayArray(k)*cosd(angleArray(j)) - pointIndex1(8)*tand(angleArray(j)));
end
if (rayArray(k)<0)
    pointIndex1(2) = yRow/2;
    pointIndex1(1) = (abs(rayArray(k))*sind(angleArray(j)+180)*tand(angleArray(j)+180) + abs(rayArray(k))*cosd(angleArray(j)+180) - pointIndex1(2)*tand(angleArray(j)+180));
    pointIndex1(3) = -xRow/2;
    pointIndex1(4) = (-cotd(angleArray(j)+180)*pointIndex1(3) + abs(rayArray(k))*sind(angleArray(j)+180) + cotd(angleArray(j)+180)*abs(rayArray(k))*cosd(angleArray(j)+180));
    pointIndex1(5) = xRow/2;
    pointIndex1(6) = (-cotd(angleArray(j)+180)*pointIndex1(5) + abs(rayArray(k))*sind(angleArray(j)+180) + cotd(angleArray(j)+180)*abs(rayArray(k))*cosd(angleArray(j)+180));
    pointIndex1(8) = -yRow/2;
    pointIndex1(7) = (abs(rayArray(k))*sind(angleArray(j)+180)*tand(angleArray(j)+180) + abs(rayArray(k))*cosd(angleArray(j)+180) - pointIndex1(8)*tand(angleArray(j)+180));
end
correctPoints = isfinite(pointIndex1);
number = 0;
for a = 1:7
    if (correctPoints(a) == 1 && correctPoints(a+1)==1)
        number = number + 1;
    end
end
pointIndex2 = zeros(1, number);
positionPoint = 1;
for x = 1:2:8
    if (isfinite(pointIndex1(x)) == 1 && isfinite(pointIndex1(x+1)))
        pointIndex2(positionPoint) = pointIndex1(x);
        pointIndex2(positionPoint+1) = pointIndex1(x+1);
        positionPoint = positionPoint + 2;
    end
end
positionPoint = 1;
pointIndex3 = zeros(1,4);
for pt = 1:2:length(pointIndex2)
    if ((pointIndex2(pt) >= pointIndex1(3) && (pointIndex2(pt) <= pointIndex1(5) && (pointIndex2(pt+1) <= pointIndex1(2) && pointIndex2(pt+1) >= pointIndex1(8)))))
        pointIndex3(positionPoint) = pointIndex2(pt);
        pointIndex3(positionPoint+1) = pointIndex2(pt+1);
        positionPoint = positionPoint + 2;
    end
end
pointIndex4 = zeros(1,4);
amount = 0;
pointIndex4(1)=pointIndex3(1);
pointIndex4(2)=pointIndex3(2);
for pos = 3:2:length(pointIndex3)-1
    if (pointIndex4(1)~= pointIndex3(pos))
        pointIndex4(3) = pointIndex3(pos) ;
        pointIndex4(4) = pointIndex3(pos+1) ;
        amount = amount + 1;
    end
end
if (amount == 0)
    pointIndex4(3) = pointIndex3(3) ;
    pointIndex4(4) = pointIndex3(4) ;
end
end
%}