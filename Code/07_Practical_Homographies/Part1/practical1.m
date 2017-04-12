function practical1

%The aim of practical 1 is to calculate the homography that best maps two
%sets of points to one another.  We will (eventually) use this for creating
%panoramas, and for calculating the 3d pose of planes.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%close all open figures
close all;

%set of two dimensional Cartesian points
pts1Cart = [  240.5000   16.8351   33.5890  164.2696  149.1911;...
              248.8770  193.5890   251.3901 168.4581  228.7723];

%turn points to homogeneous representation
pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))]
      
%define a homography
H = [0.6 0.7 -100; 1.0 0.6 50; 0.001 0.002 1.0]

%apply homography to points
pts2Hom = H*pts1Hom

%convert back to Cartesian
%repmat: repeat mat alone different axis
pts2Cart = pts2Hom(1:2,:)./repmat(pts2Hom(3,:),2,1)

%add a small amount of noise
noiseLevel = 4.0;
%pts2Cart = pts2Cart+noiseLevel*randn(size(pts2Cart));

%draw two set of two dimensional points
%opens figure
figure; set(gcf,'Color',[1 1 1]);
%draw lines between each pair of points
nPoint = size(pts1Cart,2)
for (cPoint = 1:nPoint)
    %plot a green line between each pair of points
    plot([pts1Cart(1,cPoint) pts2Cart(1,cPoint)],[pts1Cart(2,cPoint) pts2Cart(2,cPoint)],'g-');
    %make sure we don't replace with next point
    hold on;
end;

%draws first set of points
plot(pts1Cart(1,:),pts1Cart(2,:),'b.','MarkerSize',20);
%remove axis
set(gca,'Box','Off');

%draws second set of points
plot(pts2Cart(1,:),pts2Cart(2,:),'r.','MarkerSize',20);

%now our goal is to transform the first points so that they map to the
%second set of points

%****TO DO****: Fill in the details of this routine 
%At the moment, it just returns and identity matrix (body is below)
scale = 1;
HEst = calcBestHomography(pts1Cart, pts2Cart);

HEst = scale * HEst; 

%now we will see how well the routine works by applying the mapping and
%measuring the square  distance between the desired and actual positions

%apply homography to points
pts2EstHom = HEst*pts1Hom

%convert back to Cartesian
pts2EstCart = pts2EstHom(1:2,:)./repmat(pts2EstHom(3,:),2,1)

%calculate mean squared distance from actual points
sqDiff = mean(sum((pts2Cart-pts2EstCart).^2))
sqDiff
%draw figure with points before and after
%draw two set of two dimensional points
%opens figure
figure; set(gcf,'Color',[1 1 1]);
%draw lines between each pair of points
nPoint = size(pts1Cart,2)
for (cPoint = 1:nPoint)
    %plot a green line pairs of actual and estimated points
    plot([pts2Cart(1,cPoint) pts2EstCart(1,cPoint)],[pts2Cart(2,cPoint) pts2EstCart(2,cPoint)],'g-');
    %make sure we don't replace with next point
    hold on;
end;

%draws second set of points
plot(pts2Cart(1,:),pts2Cart(2,:),'r.','MarkerSize',20);
%remove axis
set(gca,'Box','Off');

%draws estimated positions of second set of points
plot(pts2EstCart(1,:),pts2EstCart(2,:),'m.','MarkerSize',20);

%other things **** TO DO ****
%1. Convince yourself that the homography is ambiguous up to scale (by
%multiplying it by a constant factor and showing it does the same thing).
%Can you see why this is mathematically the case?
%2. Show empirically that your homography routine can EXACTLY map any four points to any
%other four points
%3. Now move to practical 1b.

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


