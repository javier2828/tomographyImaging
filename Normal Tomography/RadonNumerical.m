% Javier Salazar - Mathematical Radon Part 1 - Squares - call: RadonNumerical(365,180, 1)
function [] = RadonNumerical(numberRays,numberAngles, numSquares) % mathematial radon inputs
    rayArray = linspace(-sqrt(2),sqrt(2),numberRays); % array for ray values
    angleArray = linspace(1,180,numberAngles); % angle storage
    grid = zeros(numberRays, numberAngles); % radon grid setup
    sqX = zeros(1,numSquares); %square array info
    sqY = zeros(1,numSquares);
    sqL = zeros(1,numSquares);
    sqD = zeros(1,numSquares);
    for i = 1:numSquares   % input square information
        fprintf('Please enter info for square #%d:\n',i);
        sqX(i) = input('Please enter x position [-1,1]:\n');
        sqY(i) = input('Please enter y position [-1,1]:\n');
        sqL(i) = input('Please enter side length [0,2]:\n');
        sqD(i) = input('Please enter density [0,1]:\n');
    end
    tic % start stopwatch
    for sq = 1:numSquares % fill for all squares
        for j = 1:numberAngles % all elements in angle array
            for k = 1:numberRays % all elements in rays
                    if (rayArray(k) >= 0) % if ray is positive do this
                        pointIndex = zeros(1,8); % array for 4 max point values for sq boundary      
                        filteredIndex = zeros(1,4); % filtered index of 2 points
                        pointIndex(2) = (sqL(sq)/2)+sqY(sq); % y postion of first point
                        pointIndex(1) = (rayArray(k)*sind(angleArray(j))*tand(angleArray(j)) + rayArray(k)*cosd(angleArray(j)) - pointIndex(2)*tand(angleArray(j))); % xposition of point 1
                        pointIndex(3) = -(sqL(sq)/2)+sqX(sq); % x position 2nd point
                        pointIndex(4) = (-cotd(angleArray(j))*pointIndex(3) + rayArray(k)*sind(angleArray(j)) + cotd(angleArray(j))*rayArray(k)*cosd(angleArray(j))); % y position point 2
                        pointIndex(5) = (sqL(sq)/2)+sqX(sq); % x pos of point 3
                        pointIndex(6) = (-cotd(angleArray(j))*pointIndex(5) + rayArray(k)*sind(angleArray(j)) + cotd(angleArray(j))*rayArray(k)*cosd(angleArray(j))); % y position of point 3
                        pointIndex(8) = -(sqL(sq)/2)+sqY(sq); % y position for lower row
                        pointIndex(7) = (rayArray(k)*sind(angleArray(j))*tand(angleArray(j)) + rayArray(k)*cosd(angleArray(j)) - pointIndex(8)*tand(angleArray(j))); % x pos for row intersection
                        correctPoints = isfinite(pointIndex); % checks which points are valid
                        number = 0; % counts valid points to form new array
                        for a = 1:7 % goes through which point pairs work
                            if (correctPoints(a) == 1 && correctPoints(a+1)==1)
                                number = number + 1;
                            end
                        end
                        preIndex = zeros(1, number); % new index based correct filtered points
                        positionPoint = 1; % pass point info to new array preindex
                        for x = 1:2:8
                            if (isfinite(pointIndex(x)) == 1 && isfinite(pointIndex(x+1)))
                                preIndex(positionPoint) = pointIndex(x);
                                preIndex(positionPoint+1) = pointIndex(x+1);
                                positionPoint = positionPoint + 2;
                            end
                        end
                        positionPoint = 1;
                        for pt = 1:2:length(preIndex)   % checks which points fall within boundary and pass info to filteredIndex
                             if ((preIndex(pt) >= pointIndex(3) && (preIndex(pt) <= pointIndex(5) && (preIndex(pt+1) <= pointIndex(2) && preIndex(pt+1) >= pointIndex(8)))))
                                filteredIndex(positionPoint) = preIndex(pt);
                                filteredIndex(positionPoint+1) = preIndex(pt+1);
                                positionPoint = positionPoint + 2;
                             end
                        end 
                        filteredIndex2 = zeros(1,4); % one more array for duplicate value removal
                        amount = 0; % keeps track of non duplicate
                        filteredIndex2(1)=filteredIndex(1); % point 1 not duplicate only other points can be
                        filteredIndex2(2)=filteredIndex(2);
                        for pos = 3:2:length(filteredIndex)-1
                            if (filteredIndex2(1)~= filteredIndex(pos)) % if point not duplicate pass on array
                                filteredIndex2(3) = filteredIndex(pos) ;
                                filteredIndex2(4) = filteredIndex(pos+1) ;
                                amount = amount + 1;
                            end
                        end
                        if (amount == 0) % if duplicate only point then pass values
                            filteredIndex2(3) = filteredIndex(3) ;
                            filteredIndex2(4) = filteredIndex(4) ;
                        end
                        pointLength = sqrt((filteredIndex2(3)-filteredIndex2(1))^2 + (filteredIndex2(4)-filteredIndex2(2))^2); % length between two points
                        pointWeight = pointLength*sqD(sq); % multtiply by density
                        grid(k,j) = grid(k,j) + pointWeight ; % fill in grid                 
                    end
                    if (rayArray(k) < 0) % same for negative ray values with abs(ray) and 180 degree radon shift property
                        pointIndex = zeros(1,8);
                        filteredIndex = zeros(1,4);
                        pointIndex(2) = (sqL(sq)/2)+sqY(sq);
                        pointIndex(1) = (abs(rayArray(k))*sind(angleArray(j)+180)*tand(angleArray(j)+180) + abs(rayArray(k))*cosd(angleArray(j)+180) - pointIndex(2)*tand(angleArray(j)+180));
                        pointIndex(3) = -(sqL(sq)/2)+sqX(sq);
                        pointIndex(4) = (-cotd(angleArray(j)+180)*pointIndex(3) + abs(rayArray(k))*sind(angleArray(j)+180) + cotd(angleArray(j)+180)*abs(rayArray(k))*cosd(angleArray(j)+180));
                        pointIndex(5) = (sqL(sq)/2)+sqX(sq);
                        pointIndex(6) = (-cotd(angleArray(j)+180)*pointIndex(5) + abs(rayArray(k))*sind(angleArray(j)+180) + cotd(angleArray(j)+180)*abs(rayArray(k))*cosd(angleArray(j)+180));
                        pointIndex(8) = -(sqL(sq)/2)+sqY(sq);
                        pointIndex(7) = (abs(rayArray(k))*sind(angleArray(j)+180)*tand(angleArray(j)+180) + abs(rayArray(k))*cosd(angleArray(j)+180) - pointIndex(8)*tand(angleArray(j)+180));
                        correctPoints = isfinite(pointIndex);
                        number = 0;
                        for a = 1:7
                            if (correctPoints(a) == 1 && correctPoints(a+1)==1)
                                number = number + 1;
                            end
                        end
                        preIndex = zeros(1, number);
                        positionPoint = 1;
                        for x = 1:2:8
                            if (isfinite(pointIndex(x)) == 1 && isfinite(pointIndex(x+1)))
                                preIndex(positionPoint) = pointIndex(x);
                                preIndex(positionPoint+1) = pointIndex(x+1);
                                positionPoint = positionPoint + 2;
                            end
                        end
                        positionPoint = 1;
                        for pt = 1:2:length(preIndex)
                             if (((preIndex(pt) >= pointIndex(3) && (preIndex(pt) <= pointIndex(5))) && ((preIndex(pt+1) <= pointIndex(2)) && (preIndex(pt+1) >= pointIndex(8)))))
                                filteredIndex(positionPoint) = preIndex(pt);
                                filteredIndex(positionPoint+1) = preIndex(pt+1);
                                positionPoint = positionPoint + 2;
                             end
                        end
                        filteredIndex2 = zeros(1,4);
                        amount = 0;
                        filteredIndex2(1)=filteredIndex(1);
                        filteredIndex2(2)=filteredIndex(2);
                        for pos = 3:2:length(filteredIndex)-1
                            if (filteredIndex2(1)~= filteredIndex(pos))
                                filteredIndex2(3) = filteredIndex(pos) ;
                                filteredIndex2(4) = filteredIndex(pos+1) ;
                                amount = amount + 1;
                            end
                        end
                        if (amount == 0)
                            filteredIndex2(3) = filteredIndex(3) ;
                            filteredIndex2(4) = filteredIndex(4) ;
                        end
                        pointLength = sqrt((filteredIndex2(3)-filteredIndex2(1))^2 + (filteredIndex2(4)-filteredIndex2(2))^2);
                        pointWeight = pointLength*sqD(sq);
                        grid(k,j) = grid(k,j) + pointWeight ;
                    end    
            end
        end
    end
    figure('Name','Radon ALG-1 Function','NumberTitle','off') % starts new figure for radon grid
    imagesc(grid, 'XData', [1 180], 'YData', []); % scale grid and display range of data
    axis image % keeps image proportionality
    title(['Radon ALG-1 Function ',num2str(size(grid,2)),' X ',num2str(size(grid,1))]); % title showing grid size
    toc % end stopwatch
end