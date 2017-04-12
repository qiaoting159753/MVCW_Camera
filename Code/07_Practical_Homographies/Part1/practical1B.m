function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

%close all open figures
close all;

%load in the required data
load('PracticalData','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3

%show images and points
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
plot(pts1(1,:),pts1(2,:),'r.'); 
plot(pts1b(1,:),pts1b(2,:),'m.');

figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
plot(pts2(1,:),pts2(2,:),'r.'); 

figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
plot(pts3(1,:),pts3(2,:),'m.'); 

%****TO DO**** 
%calculate homography from pts1 to pts2
H12 = calcBestHomography(pts1,pts2);

%****TO DO**** 
for row=1:size(im1,1)
    for col=1:size(im1,2)
        %transform this pixel position with your homography to find where it 
        %is in the coordinates of image 2
        transfered = H12 * [col;row;1];
        transfered = round(transfered(1:2,:)./transfered(3,:));

        %if it the transformed position is within the boundary of image 2 then 
        %copy pixel colour from image 2 pixel to current position in image 1 
        %draw new image1 (use drawnow to force it to draw)
        if (transfered(1)<size(im2,2) && transfered(1)>0 && transfered(2)<size(im2,1)&& transfered(2)>0)
            im1(row,col,:)=im2(transfered(2),transfered(1),:);
        end
    end
end;

%****TO DO****
%repeat the above process mapping image 3 to image 1.
H13 = calcBestHomography(pts1b,pts3);

%****TO DO**** 
for row=1:size(im1,1)
    for col=1:size(im1,2)
        %transform this pixel position with your homography to find where it 
        %is in the coordinates of image 2
        transfered = H13 * [col;row;1];
        transfered = round(transfered(1:2,:)./transfered(3,:));

        %if it the transformed position is within the boundary of image 2 then 
        %copy pixel colour from image 2 pixel to current position in image 1 
        %draw new image1 (use drawnow to force it to draw)
        if (transfered(1)<size(im3,2) && transfered(1)>0 && transfered(2)<size(im3,1)&& transfered(2)>0)
            im1(row,col,:)=im3(transfered(2),transfered(1),:);
        end
    end
end;

figure;imshow(uint8(im1));
%==========================================================================
function H = calcBestHomography(pts1Cart, pts2Cart)
%**** TO DO ****;
%first turn points to homogeneous
pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))];
pts2Hom = [pts2Cart; ones(1,size(pts1Cart,2))];

%xi' = H * xi
%pts2Hom = H * pts1Home
%then construct A matrix which should be (10 x 9) in size
%solve Ah = 0 by calling
A = zeros(10,9);
t_0 = zeros(1,3);

i = 1;
for index = 1:2:9
    A(index,1:3)=t_0;
    A(index+1,4:6)=[0,0,0];
    A(index+1,1:3)=[pts1Hom(1,i),pts1Hom(2,i),1];
    A(index,4:6)=[-pts1Hom(1,i),-pts1Hom(2,i),-1];
    A(index,7:9)=[pts2Hom(2,i)*pts1Hom(1,i),pts2Hom(2,i)*pts1Hom(2,i),pts2Hom(2,i)];
    A(index+1,7:9)=[-pts2Hom(1,i)*pts1Hom(1,i),-pts2Hom(1,i)*pts1Hom(2,i),-pts2Hom(1,i)];
    i = i + 1;
end



h = solveAXEqualsZero(A);
H = reshape(h,3,3)';

%reshape h into the matrix H

%Beware - when you reshape the (9x1) vector x to the (3x3) shape of a homography, you must make
%sure that it is reshaped with the values going first into the rows.  This
%is not the way that the matlab command reshape works - it goes columns
%first.  In order to resolve this, you can reshape and then take the
%transpose


%==========================================================================
function x = solveAXEqualsZero(A);
    [U,S,V]=svd(A);
    v_size = size(V);
    x = V(:,v_size(2));


