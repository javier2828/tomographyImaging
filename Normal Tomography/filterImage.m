%Javier Salazar - Filter algorithm
%call: filterImage(radonGrid, 1)
%filterOption:  1-ramp 2-cosine 3-gaussian 4-sinc 5-exponantial
function [filteredGrid] = filterImage(radonGrid, filterOption) % filter the radon image.
tic %start stopwatch
[yRow, xRow] = size(radonGrid); % get dimensions of radon image
rayArray = linspace(-yRow/2,yRow/2,yRow); % create array of t values used
Fs = abs(2*(1/(rayArray(2)-rayArray(1)))); % use ray line spacing to determine sampling frequency
%freqArray = (-Fs/2:Fs/(yRow-1):Fs/2); % ignore
freqArray = linspace(-Fs/2,Fs/2,yRow); % create frequency domain values from sampling frequency
fourierGrid = fft(radonGrid); % convert radon to frequency domain
fourierGrid = fftshift(fourierGrid,1); % shift zero frequency to center
if (filterOption == 1) % ram-lak filter
    filter = abs(freqArray');
end
if (filterOption == 2) % hann filter
    filter = abs(freqArray').*cos((2/Fs)*acos(0)*freqArray');
end
if (filterOption == 3) % gaussian filter
    filter = abs(freqArray').* exp(((4*log(0.01))/(Fs^2))*(freqArray').^2);
end
if (filterOption == 4) % sinc filter
    filter = abs(freqArray').* ((sin(2*freqArray'))./(2*freqArray'));
end
if (filterOption == 5) % exponential filter
    filter = abs(freqArray').* exp(-3*abs(freqArray'));
end
if (filterOption == 6) % exponential filter
    filter = (freqArray').^2;
end
filter = filter ./ max(filter); % normalize filter
freqGrid = repmat(filter,1,xRow); % convert filter array values to matrix since values are repeated
fourierGrid = fourierGrid.*freqGrid; % multiply grid by filter grid
fourierGrid = ifftshift(fourierGrid,1); % shift zero frequency back to edge before converting back
filteredGrid = ifft(fourierGrid,'symmetric'); % go back to spatial domain. symmetric needed to force real values due to conjugate symmetry.
figure('Name','Filter Design','NumberTitle','off') % open new figure to show filter used
plot(freqArray,filter); % plot filter values versus frequency domain
title('Filter Design'); % simple title
ylabel('Magnitude'); % label axis
xlabel('Frequency');
toc % end stopwatch
end