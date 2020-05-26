% Javier Salazar - Discrete Radon Transform Part 2 - call: imageRadon('phantom.jpg',180,1,1)
% mode = 0 - ALG-2 ONLY! mode = 1 - ALG-2 & MATLAB RADON - TWO FIGURES!
function [finalGrid] = newRadon(imageMatrix, numberAngles, mode) % new radon that works on rectangular images
tic % start stopwatch
global imageGrid; % inialize some global variables
imageGrid = imageMatrix;
global xRow;
global yRow;
global rayArray;
global angleArray;
global radonGrid;
[yRow, xRow] = size(imageGrid); % dimensions from phantom
diagPixels = ceil(sqrt(xRow^2+yRow^2))+4; % how many t values we will do +4
angleArray = linspace(1,180,numberAngles); % angles that will define the resolution. 180 works best
rayArray = linspace(-diagPixels/2,diagPixels/2,diagPixels); % ray distances
radonGrid = zeros(diagPixels,numberAngles); % start grid
for j = 1:numberAngles
    for k = 1:diagPixels
       % if (angleArray(j)==90 && k==109)
            [columnPoints,rowPoints] = getPoints(j,k); % will get points based on location of pixels
            [xPoints, yPoints] = pointFilter(columnPoints,rowPoints); % additional filtration before filling grid
            rValue = fillGrid(xPoints, yPoints); % will return grid value based on point grid values
            radonGrid(k,j) = rValue; % value assigned
        %end
    end
end
radonGrid2 = radonGrid; % second grid to fill from 180-360
radonGrid = flipud(radonGrid); % must be flipped
finalGrid = flipud(horzcat(radonGrid,radonGrid2)); % joins both grids to one
imwrite(mat2gray(finalGrid),'RadonALG2Phantom.jpg'); % write grid to file converted to greyscale
figure('Name','RADON ALG-2 IMAGE','NumberTitle','off'); % opens new figure
iptsetpref('ImshowAxesVisible','on'); % show axis numbers
imagesc(finalGrid, 'XData', [1 360], 'YData', [-diagPixels/2 diagPixels/2]); % ranges of exis and scale image to display
axis image; % keep proportionality
title(['RADON ALG-2 FUNCTION: ',num2str(numberAngles*2),' X ',num2str(diagPixels)]); % title the figure
xlabel('Angle Values'); % label axis
ylabel('Ray Values');
h = colorbar; % show density bar
ylabel(h, 'Density');
if (mode == 1) % same as above but using MATLAB RADON FUNCTION - TWO FIGURES
    radonMatrix = radon(imageGrid, angleArray);
    radonMatrix2 = flipud(radonMatrix);
    radonMatrix = horzcat(radonMatrix,radonMatrix2);
    imwrite(mat2gray(radonMatrix),'MATLABFUNCTION-Phantom.jpg');
    figure('Name','MATLAB FUNCTION IMAGE','NumberTitle','off');
    iptsetpref('ImshowAxesVisible','on');
    imagesc(radonMatrix, 'XData', [1 360], 'YData', [-diagPixels/2 diagPixels/2]);
    axis image;
    title(['RADON MATLAB FUNCTION: ',num2str(length(radonMatrix(1,:))),' X ',num2str(length(radonMatrix(:,1)))]);
    xlabel('Angle Values');
    ylabel('Ray Values');
    h = colorbar;
    ylabel(h, 'Density');
end
toc % stop stopwatch
end

function [radonValue] = fillGrid(xPoints, yPoints) % figure out radon value at point
global imageGrid; % intialize some global variables
global xRow;
global yRow;
pointLength = length(xPoints); % how many points we have
radonValue = 0;
columns = -xRow/2:xRow/2; % grid intersections 
rows = yRow/2:-1:-yRow/2;
imageDensity = zeros(1,pointLength-1); % pixel location values
for i= 1:pointLength-1
    lengthLine = sqrt( (yPoints(i+1)-yPoints(i))^2 + (xPoints(i+1)-xPoints(i))^2); % length of line intersection
    xValue = floor((xPoints(i) + (xPoints(i+1)-xPoints(i))/2)); % figure out x and y location of pixel 
    yValue = ceil((yPoints(i) + (yPoints(i+1)-yPoints(i))/2));
    columnValue = find(columns == xValue); % transform value to image grid values since 1,1 is upperleft corner
    rowValue = find(rows == yValue);
    imageDensity(i) = imageGrid(rowValue,columnValue); % get density value of pixel
    densityLine = lengthLine*imageDensity(i); % multiply by length of line
    radonValue = radonValue + densityLine; % add for all points
end
end

function [sortedX, sortedY] = pointFilter(columnPoints,rowPoints) % bring in points
columnPoints = columnPoints(1:find(columnPoints,1,'last')); % remove any trailing zeros left over
rowPoints = rowPoints(1:find(rowPoints,1,'last'));
rowLength = length(rowPoints); % get length of points
columnLength = length(columnPoints);
columnPoints2 = columnPoints;
rowPoints2 = rowPoints;
if (mod(columnLength,2)==1) % in case a trailing zero was actually a value from points alg. then pass a zero at the end
    columnPoints2 = zeros(1,columnLength+1);
    columnPoints2(1:columnLength) = columnPoints(1:columnLength);
    columnLength = columnLength + 1;
end
if (mod(rowLength,2)==1) % same as above
    rowPoints2 = zeros(1,rowLength+1);
    rowPoints2(1:rowLength) = rowPoints(1:rowLength);
    rowLength = rowLength + 1;
end
columnPoints3 = zeros(1,columnLength); % new column points to remove duplicate points
xPos = 1;
number = 0;
totalNumber = 0;
for x = 1:2:columnLength % go through all of columnpoints and compare against row points
    for k = 1:2:rowLength
        if (columnPoints2(x) == rowPoints2(k)) % keeps track of duplicates
           number = number + 2;
        end
    end
    if (number == 0) % if its a non duplicate value then it will be passed along new array
            columnPoints3(xPos) = columnPoints2(x);
            columnPoints3(xPos+1) = columnPoints2(x+1);
            xPos = xPos + 2;
    end
    totalNumber = totalNumber + number; % some counters for length
    number = 0;
end
points = zeros(1,rowLength+columnLength-totalNumber); % merge row and column points to one array
for i = 1:rowLength
    points(i) = rowPoints2(i);
end
for i = 1:columnLength-totalNumber
   points(rowLength+i) = columnPoints3(i); 
end
xPoints = zeros(1,(rowLength+columnLength-totalNumber)/2); % start new arrays to keep x values in one and y values in the other
yPoints = zeros(1,(rowLength+columnLength-totalNumber)/2);
xPos = 1;
yPos = 1;
for i = 1:2:rowLength+columnLength-totalNumber % pass values
    xPoints(xPos) = points(i);
    xPos = xPos + 1;
end
for i = 2:2:rowLength+columnLength-totalNumber
    yPoints(yPos) = points(i);
    yPos = yPos + 1;
end
[sortedX,sortIndex] = sort(xPoints); % sort x points from smallest to biggest
sortedY = yPoints(sortIndex); % sort y points according to the index from the x values
end

function [xIndex3,yIndex3] = getPoints(j,k) % get column and row points with some filters
global xRow; % initialize a few variables
global yRow;
global rayArray;
global angleArray;
columns = -xRow/2:0.5:xRow/2; % grid intersections
rows = -yRow/2:0.5:yRow/2;
yIndex1 = zeros(1,yRow*2+2); % max values possible
xIndex1 = zeros(1,xRow*2+2);
if(rayArray(k) >= 0) % use equations for nonnegative line values
    for i=1:2:yRow*2+1 % go trough and fill x points from row information
        yIndex1(i) = (rayArray(k)*sind(angleArray(j))*tand(angleArray(j)) + rayArray(k)*cosd(angleArray(j)) - rows(i)*tand(angleArray(j)));
        yIndex1(i+1) = rows(i);
    end
    for x=2:2:xRow*2+2 % go through and fill y points given column information
        xIndex1(x) = (-cotd(angleArray(j))*columns(x-1) + rayArray(k)*sind(angleArray(j)) + cotd(angleArray(j))*rayArray(k)*cosd(angleArray(j)));
        xIndex1(x-1) = columns(x-1);
    end
end
if (rayArray(k) < 0) % use shift property for negative t values
    for i=1:2:yRow*2+1
        yIndex1(i) = (abs(rayArray(k))*sind(angleArray(j)+180)*tand(angleArray(j)+180) + abs(rayArray(k))*cosd(angleArray(j)+180) - rows(i)*tand(angleArray(j)+180));
        yIndex1(i+1) = rows(i);
    end
    for x=2:2:xRow*2+2
        xIndex1(x) = (-cotd(angleArray(j)+180)*columns(x-1) + abs(rayArray(k))*sind(angleArray(j)+180) + cotd(angleArray(j)+180)*abs(rayArray(k))*cosd(angleArray(j)+180));
        xIndex1(x-1) = columns(x-1);
    end    
end
correctxPoints = isfinite(xIndex1); % shows where we have NaN values
correctyPoints = isfinite(yIndex1);
numberX = 0; % simple counters
numberY = 0;
for a = 1:2:xRow*2+1
   if (correctxPoints(a) == 1 && correctxPoints(a+1)==1) % count how many good values we have. NaN values discarded
       numberX = numberX + 2;
   end
end
for a = 1:2:yRow*2+1 % same for yIndex
   if (correctyPoints(a) == 1 && correctyPoints(a+1)==1)
       numberY = numberY + 2;
   end
end
xIndex2 = zeros(1, numberX); % start new index to pass non NaN values
yIndex2 = zeros(1, numberY);
positionxPoint = 1;
positionyPoint = 1;
for x = 1:2:numberX
   if (isfinite(xIndex1(x)) == 1 && isfinite(xIndex1(x+1))==1) % pass good non-NaN values from old index to new one
         xIndex2(positionxPoint) = xIndex1(x);
         xIndex2(positionxPoint+1) = xIndex1(x+1);
         positionxPoint = positionxPoint + 2;
   end
end
for x = 1:2:numberY
   if (isfinite(yIndex1(x)) == 1 && isfinite(yIndex1(x+1))) % do the same from yIndex
         yIndex2(positionyPoint) = yIndex1(x);
         yIndex2(positionyPoint+1) = yIndex1(x+1);
         positionyPoint = positionyPoint + 2;
   end
end
xValues = positionxPoint-1;
yValues = positionyPoint-1;
xIndex3 = zeros(1, xValues); % start new index to remove invalid points that are out of bound
yIndex3 = zeros(1, yValues);
positionxPoint = 1;
positionyPoint = 1;
for pt = 1:2:xValues-1   
     if (xIndex2(pt+1) <= yRow/2 && xIndex2(pt+1) >= -yRow/2 && xIndex2(pt) <= xRow/2 && xIndex2(pt) >= -xRow/2) % make sure its within the box range of image
         xIndex3(positionxPoint) = xIndex2(pt);
         xIndex3(positionxPoint+1) = xIndex2(pt+1);
         positionxPoint = positionxPoint + 2;
     end
end
for pt = 1:2:yValues-1 % same but for row points
     if (yIndex2(pt) <= xRow/2 && yIndex2(pt) >= -xRow/2 && yIndex2(pt+1) <= yRow/2 && yIndex2(pt+1) >= -yRow/2)
         yIndex3(positionyPoint) = yIndex2(pt);
         yIndex3(positionyPoint+1) = yIndex2(pt+1);
         positionyPoint = positionyPoint + 2;
     end
end
end