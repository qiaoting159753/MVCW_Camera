function r=practical2

%This project explores the geometry of a single camera. The aim is to take several points on
%a plane, and predict where they will appear in the camera image. Based on these observed
%points, we will then try to re-estimate the Euclidean transformation relating the plane and
%the camera. In practical 2b we will use this code to draw a wireframe cube
%on an augmented reality marker.   You should use this
%template for your code and fill in the missing sections marked "TO DO"


%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];
 
%We will assume an object co-ordinate system with the Z-axis pointing upwards and the
%origin in the centre of the plane. There are four known points on the plane, with coordinates
%(in mm):

XCart = [-100 -100  100  100 0 ;...
         -100  100  100 -100 0;...
          0    0    0    0   0 ];

%We will assume that the correct transformation from the plane co-ordinate system to the
%camera co-ordinate system (extrinsic matrix) is:

T = [ 0.9851  -0.0492  0.1619  46.00;...
     -0.1623  -0.5520  0.8181  70.00;...
      0.0490  -0.8324 -0.5518  500.89;...
      0        0       0       1]
  
% TO DO  Use the general pin-hole projective camera model discussed in the lectures to estimate 
%where the four points on the plane will appear in the image.  Fill in the
%details of the function "projectiveCamera" - body of function appears below

xImCart = projectiveCamera(K,T,XCart);

% TO DO Add noise to the pixel positions to simulate having to find these points in a noisy
%image. Store the results back in xImCart.  
%The noise should have standard deviation of one pixel in each direction.
xImCart = xImCart+1*randn(size(xImCart));

%Now we will take the image points and the known positions on the card and try to
%estimate the extrinsic matrix using the algorithm discussed in the lecture. 
%Fill in the details of the function "estimate plane pose" - body of function appears
%below

TEst = estimatePlanePose(xImCart,XCart,K)

%if you have got this correct, it should resemble T above.

%==========================================================================
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
alpha = sum(sum(R_top(:,1:2)./H(:,1:2)))/6;
t = (1/alpha).*H(:,4);

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
    A(i*2,9:12) = -Y(i)*pts1Hom(:,i);
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

