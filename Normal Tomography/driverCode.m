% Javier Salazar - Driver Main Function (Mother Function)
function [backGrid] = driverCode() % call: driverCode(); 
imageGrid = Phantom(128,128,0,0,1,0); % fill in inputs here for phantom
radonGrid = newRadon(imageGrid, 180, 1); % fill in inputs for radon
filteredGrid = filterImage(radonGrid, 1); % filter and enhance image
backGrid = backProjection(filteredGrid, 'cubic', 1, imageGrid); % fill in inputs. unfiltered back projection.
end