% Javier Salazar - Phantom Creator - 05/18/2017 - Call Phantom(xRow,yRow,mode,numX) on command prompt;
% xRow = number of horizontal pixels, enter even number.
% yRow = number of vertical pixels, enter even number.
% Mode = 0 - overlay last shape to grid, 1 - merge density values
% numX = number of specific shapes
function [grid2] = Phantom(xRow,yRow, mode, numCircle, numRect, numTri) % main function that draws phantom
global grid; % global so that other funtions can write to grid
grid = zeros(yRow,xRow); % the grid for the phantom
circleX = zeros(1,numCircle); % Array to store x values of circes
circleY = zeros(1,numCircle); % same as above for y values
circleR = zeros(1,numCircle); % circle radici
circleD = zeros(1,numCircle); % density values
rectX = zeros(1,numRect); % same for rectangles X
rectY = zeros(1,numRect); % Y positions
rectL = zeros(1,numRect); % length values
rectH = zeros(1,numRect); % height values
rectD = zeros(1,numRect); %density values
triX = zeros(1,numTri); % triangle  position
triY = zeros(1,numTri); % triangle y positions
triL = zeros(1,numTri); % length of one side of triangle
triD = zeros(1,numTri); % density values
for i = 1:numCircle   % fills array information based on user input
    fprintf('Please enter info for circle #%d:\n',i);
    circleX(i) = input('Please enter x position [-1,1]:\n'); % input info
    circleY(i) = input('Please enter y position [-1,1]:\n');
    circleR(i) = input('Please enter radius [0,1]:\n');
    circleD(i) = input('Please enter density [0,1]:\n');
    drawCircle(xRow/2,yRow/2,circleX(i),circleY(i),circleR(i),circleD(i), mode); % draws after user gives info
end
for i = 1:numRect % same as above
    fprintf('Please enter info for rectangle #%d:\n',i);
    rectX(i) = input('Please enter x position [-1,1]:\n'); % user input
    rectY(i) = input('Please enter y position [-1,1]:\n');
    rectL(i) = input('Please enter length [0,2]:\n');
    rectH(i) = input('Please enter height [0,2]:\n');
    rectD(i) = input('Please enter density [0,1]:\n');
    drawRect(xRow/2,yRow/2,rectX(i),rectY(i),rectL(i), rectH(i), rectD(i), mode); %draws after user gives info
end
for i = 1:numTri  % same as above
    fprintf('Please enter info for eq. triangle #%d:\n',i);
    triX(i) = input('Please enter x position [-1,1]:\n'); % user input
    triY(i) = input('Please enter y position [-1,1]:\n');
    triL(i) = input('Please enter length [0,2]:\n');
    triD(i) = input('Please enter density [0,1]:\n');
    drawTri(xRow/2,yRow/2,triX(i),triY(i),triL(i), triD(i), mode); % draws after user gives info
end
tic % start stopwatch
grid2 = grid ./ max(grid(:)); % value that the function returns normalized
imwrite(grid2,'phantom.jpg'); % write grid to image for future viewing
figure('Name','Phantom','NumberTitle','off') % starts new figure for phantom
imagesc(grid2, 'XData', [-1 1], 'YData', [-1 1]);
axis image % keeps image proportionality
title(['Phantom: ',num2str(xRow),' X ',num2str(yRow)]); % title information for the phantom
xlabel('X Position'); % label axis
ylabel('Y Position');
h = colorbar; % density bar
ylabel(h, 'Density'); % label density
toc % end stopwacth
end

function [] = drawTri(centerX, centerY, xPos, yPos, lengthT, tDen, mode) % necessary triangle info
global grid;
tLen = round((lengthT*centerX)/2); % half the length of a side
posX = round(centerX*(1+xPos)); % converts coordinate [0-1] to pizel grid position
posY = round(centerY*(1-yPos));
pointAX = posX-tLen; % must find position of 3 triangle points
pointAY = posY+round(tLen*2*0.288675); % y position of point 1
pointBX = posX+tLen; % barycentric triangle coordinates to check if points are within triangle
pointCY = posY-round(tLen*2*0.57735); % more point info
for k = pointCY:pointAY % searches from top y point down to bottom y point
    for j = pointAX:pointBX % starts from left most x point to right point
        as_x = j-pointAX; %relative point info
        as_y = k-pointAY; % relative y point info
        s_ab = (pointBX-pointAX)*as_y-(pointAY-pointAY)*as_x > 0; % condition whether within triangle
        if (((posX-pointAX)*as_y-(pointCY-pointAY)*as_x > 0 ~= s_ab) && ((posX-pointBX)*(k-pointAY)-(pointCY-pointAY)*(j-pointBX) > 0 == s_ab)) %barycentric coordinate formulas
            if (j>0) && (j<=centerX*2) && (k>0) && (k<=centerY*2) && (mode==0) % only writes if points are within grid, crops other points off
                grid(k,j) = tDen; % if mode 0 then just write density
            end
            if (j>0) && (j<=centerX*2) && (k>0) && (k<=centerY*2) && (mode==1)
                grid(k,j) = grid(k,j) + tDen; % otherwise superimpose info
            end
        end
    end
end
end

function [] = drawRect(centerX, centerY, xPos, yPos, lengthR, heightR, rDen, mode) % draw rectangle function
global grid;
rLen = round((lengthR*centerX)/2); % convert to grid coordinates
rHei = round((heightR*centerY)/2); % same as above
posX = round(centerX*(1+xPos)); % change position of points
posY = round(centerY*(1-yPos)); % grid point conversion
for j = posX-rLen:posX+rLen % start from one end point to other
    for k = posY-rHei:posY+rHei % start frin top to bottom
        if (j>0) && (j<=centerX*2) && (k>0) && (k<=centerY*2) && (mode==0) % if point within grid then fill info
            grid(k,j)= rDen; % mode 0 means only density
        end
        if (j>0) && (j<=centerX*2) && (k>0) && (k<=centerY*2) && (mode==1)
            grid(k,j)= grid(k,j)+rDen; % otherwise superimpose
        end
    end
end
end

function [] = drawCircle(centerX,centerY,xPos, yPos, cRad, cDen, mode) % draw circles
global grid;
posX = round(centerX*(1+xPos)); % convert x and y position
posY = round(centerY*(1-yPos));
radLength = round(min(centerY,centerX)*cRad);%convert to pixel infp
for j = posX-radLength:posX+radLength % goes through square x position
    for k = posY-radLength:posY+radLength % goes through square with circle inscrimbed
        posLength = sqrt(((posX-j)^2)+(posY-k)^2); % equation for circle
        if (posLength <= radLength) && (j>0) && (j<=centerX*2) && (k>0) && (k<=centerY*2) && (mode==0)
            grid(k,j)= cDen; % if within limits and within radius then draw values
        end
        if (posLength <= radLength) && (j>0) && (j<=centerX*2) && (k>0) && (k<=centerY*2) && (mode==1)
            grid(k,j)= grid(k,j)+cDen; % if mode 1 then superimpose
        end
    end
end
end