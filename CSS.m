close all, clear all, clc;

rgbim = imread('lesion1.jpg'); % Read in image
grayim = double(rgb2gray(rgbim))/255; % Convert to grayscale

css = []; % Initialize CSS image

limit = 20; % Set sigma limit for non malignant lesions

figure, hold
xlabel('t');
ylabel('sigma');

structurally_malignant = false;

for sigma = 0:29 % Set sigma level for gaussian blurring

    edgeim = edge(grayim,'canny',[0.4 0.6], (4+sigma)*sqrt(2)); % Edge detection with base sigma value of 4*sqrt(2)
    edgeim = filledgegaps(edgeim,7); % Fill gaps up to 7 pixels wide
    
    %subplot(5,6,sigma+1), subimage(edgeim), axis off
    %imshow(edgeim)
    
    % Create pixel wide gap at 12 o'clock for consistent edge
    for i = 1:fix(size(edgeim,1)/2)
        edgeim(i,fix(size(edgeim,2)/2)) = 0;
    end
    
    [edgelist, edgeim] = edgelink(edgeim, 1); % Convert to list of points
    edgelist = edgelist{1,1}; % Extract list from cell structure

    if sigma == 0
        max_length = double(size(edgelist,1)); % Set max length to length of edge at sigma 0
    end
    
    K = []; % Initialize curvature array
    
    pixel_jump = 20.0; % Set amount of pixels to look behind/ahead for calculating curvature
    
    % Calculate discrete curvature for every point in the edge list
    for i = pixel_jump+1:size(edgelist,1)-pixel_jump        
        x_pre = double(edgelist(i-pixel_jump,2));
        x = double(edgelist(i,2));
        x_post = double(edgelist(i+pixel_jump,2));
        
        y_pre = double(edgelist(i-pixel_jump,1));
        y = double(edgelist(i,1));
        y_post = double(edgelist(i+pixel_jump,1));
        
        dx_dt = (x_post-x_pre)/(2*pixel_jump);
        dy_dt = (y_post-y_pre)/(2*pixel_jump);
        d2x_dt2 = (x_post-(2*x)+x_pre)/(pixel_jump^2);
        d2y_dt2 = (y_post-(2*y)+y_pre)/(pixel_jump^2);
                       
        K(i-pixel_jump) = (dx_dt*d2y_dt2 - d2x_dt2*dy_dt)/((dx_dt^2 + dy_dt^2)^(3/2));
    end
    
    %plot(1:length(K),K)
    
    curr_length = double(length(K)); % Set the length of the current sigma level
    
    %figure, imshow(edgeim), hold,
    for i = 1:length(K)-1 % For every point in the list
        %plot(edgelist(i,2),edgelist(i,1),'x','LineWidth',1,'Color','r'); % Highlight pixel
        if ((sign(K(i)) == 1 && sign(K(i+1)) == -1) || (sign(K(i)) == -1 && sign(K(i+1)) == 1) || K(i) == 0) && abs(K(i)-K(i+1)) > 0.001 % If sign of next is opposite
            css(sigma+1,i) = 1;
            plot(fix(i*(max_length/curr_length)),sigma,'.','LineWidth',1,'Color','b');
            
            if sigma > limit
                structurally_malignant = true;
            end
                        
        else
            css(sigma+1,i) = 0;   
        end
    end
    
end

hline = refline([0 limit]); % Plot line to denote limit
hline.Color = 'r';

structurally_malignant