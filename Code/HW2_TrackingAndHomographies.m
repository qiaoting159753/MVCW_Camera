function HW2_TrackingAndHomographies

LLs = HW2_Practical9c( 'll' );
LRs = HW2_Practical9c( 'lr' );
ULs = HW2_Practical9c( 'ul' );
URs = HW2_Practical9c( 'ur' );

close all;

% Load frames from the whole video into Imgs{}.
% This is really wasteful of memory, but makes subsequent rendering faster.
LoadVideoFrames

% Coordinates of the known target object (a dark square on a plane) in 3D:
XCart = [-50 -50  50  50;...
          50 -50 -50  50;...
           0   0   0   0];

% These are some approximate intrinsics for this footage.
K = [640  0    320;...
     0    512  256;
     0    0    1];

% Define 3D points of wireframe object.
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];
 
hImg = figure;
       
% ================================================
for iFrame = 1:numFrames
    xImCart = [LLs(iFrame,:)' ULs(iFrame,:)' URs(iFrame,:)' LRs(iFrame,:)'];
    xImCart = circshift( xImCart, 1);

    % To get a frame from footage 
    im = Imgs{iFrame};

    % Draw image and 2d points
    set(0,'CurrentFigure',hImg);
    set(gcf,'Color',[1 1 1]);
    imshow(im); axis off; axis image; hold on;
    plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',15);


    %TO DO Use your routine to calculate TEst the extrinsic matrix relating the
    %plane position to the camera position.
    T = estimatePlanePose(xImCart, XCart, K);



    %TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
    %through the projective camera, and drawing lines betweeen the 
    %resulting 2d image points
    XImWireFrameCart = projectiveCamera(K,T,XWireFrameCart);
    hold on;
    
    % TO DO: Draw a wire frame cube using data XWireFrameCart. You need to
    % 1) project the vertices of a 3D cube through the projective camera;
    % 2) draw lines betweeen the resulting 2d image points.
    % Note: CONDUCT YOUR CODE FOR DRAWING XWireFrameCart HERE
    num_point = size(XWireFrameCart,2)
    for i = 1:num_point
        for j = 1:num_point
            %Plot a line between two points that have distance of 100 (is one edge of the cube)
            if  ( sum(abs(XWireFrameCart(:,i) - XWireFrameCart(:,j)))==100 ) 
            plot([XImWireFrameCart(1,i) XImWireFrameCart(1,j)],[XImWireFrameCart(2,i) XImWireFrameCart(2,j)],'g-');
            hold on;
            end
        end
    end
    
    hold off;
    drawnow;
    
    
end % End of loop over all frames.
% ================================================

% TO DO: QUESTIONS TO THINK ABOUT...

% Q: Do the results look realistic?
% If not then what factors do you think might be causing this


% TO DO: your routines for computing a homography and extracting a 
% valid rotation and translation GO HERE. Tips:
%
% - you may define functions for T and H matrices respectively.
% - you may need to turn the points into homogeneous form before any other
% computation. 
% - you may need to solve a linear system in Ah = 0 form. Write your own
% routines or using the MATLAB builtin function 'svd'. 
% - you may apply the direct linear transform (DLT) algorithm to recover the
% best homography H.
% - you may explain what & why you did in the report.
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
