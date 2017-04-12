function practical2b

%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%load in image 
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
             308.9825  236.7646  255.4416  340.7335   281.5895];
         
%define 3D points of plane
XCart = [-50 -50  50  50 0 ;...
          50 -50 -50  50 0;...
           0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',10);
       
%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.

TEst = estimatePlanePose(xImCart,XCart,K);
xImCart = projectiveCamera(K,TEst,XCart);


%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];
                
XImWireFrameCart = projectiveCamera(K,TEst,XWireFrameCart);
%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points
nPoint = size(XWireFrameCart,2)
for i = 1:nPoint
    for j = 1:nPoint
        %Plot a line between two points that have distance of 100 (is one edge of the cube)
        if  ( sum(abs(XWireFrameCart(:,i) - XWireFrameCart(:,j)))==100 ) 
            plot([XImWireFrameCart(1,i) XImWireFrameCart(1,j)],[XImWireFrameCart(2,i) XImWireFrameCart(2,j)],'g-');
            hold on;
        end
    end
end;



%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?
%==========================================================================



%goal of function is to project points in XCart through projective camera
%defined by intrinsic matrix K and extrinsic matrix T.
function xImCart = projectiveCamera(K,T,XCart);

%replace this
%TO DO convert Cartesian 3d points XCart to homogeneous coordinates XHom
XHom = [XCart;ones(1,size(XCart,2))];

%TO DO apply extrinsic matrix to XHom to move to frame of reference of
%camera
XCamHom = T*XHom;

%TO DO project points into normalized camera coordinates xCamHom by (achieved by
%removing fourth row)
XCamHom(4,:) = [];

%TO DO move points to image coordinates xImHom by applying intrinsic matrix
XImgHom = K*XCamHom;

%TO DO convert points back to Cartesian coordinates xImCart
xImCart = XImgHom(1:2,:)./repmat(XImgHom(3,:),2,1);

%==========================================================================
%==========================================================================

%goal of function is to estimate pose of plane relative to camera
%(extrinsic matrix) given points in image xImCart, points in world XCart
%and intrinsic matrix K.

function T = estimatePlanePose(xImCart,XCart,K)

%replace this
%TO DO Convert Cartesian image points xImCart to homogeneous representation
xImHom = [xImCart;ones(1,size(xImCart,2))];

%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
%xCamHom
xCamHom = inv(K)*xImHom;

%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.
H = calcBestHomography(XCart, xCamHom);

%TO DO Estimate first two columns of rotation matrix R from the first two
%columns of H using the SVD
%Least square sense
R = H(:,1:2);
[U,S,V] = svd(R);
R_top = U*[1,0;0,1;0,0]*V';


%TO DO Estimate the third column of the rotation matrix by taking the cross
%product of the first two columns
R_top = [R_top cross(R_top(:,1),R_top(:,2))];

%TO DO Check that the determinant of the rotation matrix is positive - if
%not then multiply last column by -1.
r_size = size(R_top);
if(det(R_top)<=0)
    R_top(:,r_size(2)) = -1*R_top(:,r_size(2));
end

%TO DO Estimate the translation t by finding the appropriate scaling factor k
%and applying it to the third colulmn of H
alpha = sum(sum(R_top(:,1:2)./H(:,1:2)));
t = alpha .* H(:,4);
t = t/6;

%TO DO Check whether t_z is negative - if it is then multiply t by -1 and
%the first two columns of R by -1.
t_size = size(t);
t_z = t(t_size(1));
if (t_z<0)
    t = -1*t;
    R_top(:,1:2) = -1*R_top(:,1:2);
end

%assemble transformation into matrix form
T  = [R_top t;0 0 0 1];

function H = calcBestHomography(pts1Cart, pts2Cart)
%**** TO DO ****;
%first turn points to homogeneous
pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))];
pts2Hom = [pts2Cart; ones(1,size(pts1Cart,2))];
%Both are 4X5 now
%So the homography matrix is 4X4

%then construct A matrix which should be (10 x 9) in size
%solve Ah = 0 by calling

[m,n] = size(pts1Hom);
%5 data points, so the size of A is (5*2)X(3*4) 
A = zeros(n * 2,3 * m);

%For each data point
%1,2,3,4,5
Y = pts2Hom(2,:);
X = pts2Hom(1,:);

for i = 1:n
    %1,3,5,7,9
    A(i*2-1,1:4) = pts1Hom(:,i);
    A(i*2-1,9:12) = -X(i)*pts1Hom(:,i);
    A(i*2,5:8) = pts1Hom(:,i);
    A(i*2,9:12) = -Y(i)*pts1Hom(:,i)
end

h = solveAXEqualsZero(A);
H = reshape(h,4,3)';

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

