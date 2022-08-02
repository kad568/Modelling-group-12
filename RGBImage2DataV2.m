function [data] = RGBImage2DataV2 (imgName)
%% Function to convert an image from with a coloured data-line into a data set automatically ...
% in the form of [data] = RGBImage2Data (img)


%% Image Data Manipulation
img = imread([imgName '.jpg']); %reads the image, specify the file type

[R,G,B]=imsplit(img);% splits image into colour channels

tholdR = 200;% defines a threshold intensity for the Red Coloour
thold = 100;% defines a threshold itnesity for Green and Blue

mask = R>tholdR & G<thold & B<thold; %creates a colour mask to filter red


%Set non-red to white
R(~mask) = 255;
G(~mask) = 255;
B(~mask) = 255;

red = cat(3,R,G,B); %creates red image


[IndexMatRed] = rgb2ind(red,2); %Converts the red image into indicies

[IndexMatBl] = rgb2ind(img,2); %Converts the black image into indicies

[y,x] = find(IndexMatRed==0); %Finds all the Red and stores the pixel coordinates of it

    %% Initialise Variables
    timeData = [];
    tempData = [];

    %% find the orgin
    [yblack,xblack] = find(IndexMatBl==0);
    origin = [mode(xblack), mode(yblack)];

    %% find the axis
    %y-Axis
    yAxVal = diff(yblack)==1; %create a logic array of where consecutive numbers exist in the y axis
    yAxVal = double(yAxVal); %convert boolean to double
    [~, ~, ~, leny] = findseq(yAxVal); %find the length of the largest string of black pixels

    yScale = 2000 / max(leny); %divide 2000 by the length of the axis in pixels to get degrees in F per pixel
    
    %x-axis
    xAxVal = xblack(find(yblack==(origin(1,2)))); %create a logic array of x coordinates for where black pixels are on the x axis
    
    xScale = 2000 /(max(xAxVal) - origin(1,1));%divide 2000 by the length of the axis in pixels to get time in S per pixel
    %% Add all the pixel data to vectors
    n=1;
    while n<=length(x)
        
        %scaled with the origin at (50,137) and 1 pixel = (7.35 deg F, 21.28s)
        timeData = [timeData, (x(n)-origin(1,1))*xScale]; 
        tempData = [tempData, (origin(1,2)-y(n))*yScale];

        n=n+1;
    end

    %Remove duplicates caused by muti-pixel thick data lines
    
    [timeData, i] = unique(timeData);
    tempData = tempData(i);

    %set last data points to 0 gradient so extrapolation is not mislead
    
    tempData(end)= tempData(end - 1);
    
    %Create variable called data for ease of use

    data = [tempData; timeData];

    %Convert F-K for SI Unit Manipulation

    data(1,:) = (data(1,:) - 32)*(5/9) + 273.15;
    
end


