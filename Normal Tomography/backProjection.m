% Javier Salazar - Backprojection Algorithm
% call: backProjection(radonGrid,'cubic',1)
% interType - 'cublic' interpolation, 'linear', 'spline', 'gauss' etc..
% mode - 0 <-- ONLY backProjection 1 fig. mode - 1 <-- my alg. & iRadon MATLAB!
function [backGrid] =  backProjection(radonGrid, interType, mode, phantomGrid) % phantomGrid imported to apply radon and iradon
tic % start stopwatch
global angleArray; % global variables needed for getLine function
global column;
global row;
global xRow;
global yRow;
global thetaRow;
[lineRow, thetaRow] = size(radonGrid); % getting info from radon input
xRow = floor((lineRow-4)/sqrt(2)); % figuring out image dimension from radonlinerow-4
yRow = xRow; % must be square dimension (iRadon does the same thing)
angleArray = linspace(1,180,thetaRow/2); % angle array 1-180
column = -xRow/2:1:xRow/2; % column locations of image grid
row = yRow/2:-1:-yRow/2; % row locations
backGrid = zeros(xRow,yRow); % initialize grid
rayGrid = linspace(-lineRow/2,lineRow/2,lineRow); % ray values
angleGrid = linspace(1,360,thetaRow); % full angles from the radon grid  
[X, Y] = meshgrid(angleGrid,rayGrid); % needed to describe mesh grid from 
for j = 1:xRow % fill in image
    for k = 1:yRow 
            [tValues, tAngles] = getLine(j,k); % get tvalues for location along with the associated angles
            interValue = interp2(X,Y,radonGrid,tAngles,tValues,interType); % interpolate using mesh grid and radonGrid
            integralMean = mean(interValue,'omitnan'); % remove any NaN values during average
            backGrid(k,j) = integralMean;  % NaN occurs from interpolation at edges (out of range) and at certain locations where angle crosses orign. very few
    end
end
backGrid(backGrid<0)=0;
backGrid = backGrid ./ max(backGrid(:)); % normalize again to compare with phantom
imwrite(mat2gray(backGrid),'backProjection-Phantom.jpg'); % write grid to jpeg
figure('Name','BACK PROJECTION IMAGE','NumberTitle','off'); % open figure
iptsetpref('ImshowAxesVisible','on'); % show axis numbers
imagesc(backGrid, 'XData', [1 xRow], 'YData', [1 yRow]); % describe axis info and scale image
axis image; % keep proportionality
title(['BACK PROJECTION FUNCTION USING ',upper(interType),' INTERPOLATION : ',num2str(xRow),' X ',num2str(yRow)]); % title the figure
xlabel('X Position'); % label axis
ylabel('y Position');
h = colorbar; % show colorbar and density values
ylabel(h, 'Density');
if (mode == 1) % perform the same to compare with iradon with no filter and the same interpolation
    radonMatrix = radon(phantomGrid, angleArray);
    backMatrix = iradon(radonMatrix, angleArray,interType);
    imwrite(mat2gray(backMatrix),'iRADON-FUNCTION-Phantom.jpg');
    figure('Name','iRADON FUNCTION IMAGE','NumberTitle','off');
    iptsetpref('ImshowAxesVisible','on');
    imagesc(backMatrix, 'XData', [1 length(backMatrix(1,:))], 'YData', [1 length(backMatrix(:,1))]);
    axis image;
    title(['iRADON FUNCTION USING ',upper(interType),' INTERPOLATION : ',num2str(length(backMatrix(1,:))),' X ',num2str(length(backMatrix(:,1)))]);
    xlabel('X Position');
    ylabel('Y Position');
    h = colorbar;
    ylabel(h, 'Density');
end
toc % stop stopwatch
end

function [tValues, tAngles] = getLine(j,k) % input grid location and output the tvalue info
global thetaRow; % must initiliaze global variables within each function
global column; 
global row;
global angleArray;
tValues = zeros(1,thetaRow/2); % initialize the array of line info
tAngles = zeros(1,thetaRow/2);
for i = 1:thetaRow/2 % go through for all angles on point 
    tValues(i) = column(j)*cosd(angleArray(i)) + row(k)*sind(angleArray(i)) ;
    tAngles(i) = angleArray(i);
    if (tValues(i) < 0)
        tValues(i) = abs(tValues(i));
        tAngles(i) = tAngles(i) + 180;
    end
end
end