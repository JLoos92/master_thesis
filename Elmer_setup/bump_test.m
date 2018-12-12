clear all
clc
figure(5)
M = zeros(300,300); % initial matrix
%N = 1; % number of bumps
sigma = 10;% std (width) of Gauss 
maxAmplitude = 10; % maximum height
[x,y] = meshgrid(1:size(M,1),1:size(M,2));
%for k=1:N
    % random location of bumps or for grounding line
    xc = 50;
    yc = 50;
    % Gauss function
    exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma^2);
    amplitude = rand()*maxAmplitude;  
    % add Gauss to the matrix M
    M = M + amplitude*exp(-exponent);
%end
surf(M)
shading interp

%%
clear allw
%Set up Cube

%(please chosse deltaX==deltaY, otherwise confusion below)
nRows=25;
nCols=nRows;
deltaX=1.0;
deltaY=deltaX;
lengthRows=nRows*deltaX;
lengthCols=nRows*deltaY;

%Set up true Velocity Structure in Block
%(we do not worry about units here. Assume it is meters per seconds)
%(use imagesc to visualiz true and guess velocity structures)
%--------------------------------------------------------------------------
VelocityInBlock=zeros(nRows,nCols)+1.0;
[X,Y]=meshgrid(1:nRows,1:nCols);
X=X/nRows;
Y=Y/nCols;
   %This is a 2D Gauss-Bell function
   %Rotate it with theta,
   %Spread in x/y defined by sigma_x/sigma_y
   %Amplitude defined by A
   %x0/y0 offset of peak in x,y direction
   theta=0*pi/180;
   sigma_x =2;
   sigma_y=2;
   x0=0.5;
   y0=0.5;
   A=100.0;
   a = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
   b = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
   c = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
   Z = A*exp( - (a*(X-x0).^2 - 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

   
   %This reshapes the square block into one model vector m for Gm=d
   %It also transforms velocity to slowness (deltaX==deltaY)
   %This needs to be done coherently, think about it what happens here.
  % m0 = reshape(deltaX./InitVelocityInBlock',nRows*nCols,1);
  % mTrue = reshape(deltaX./VelocityInBlock',nRows*nCols,1);

   %figure(1)
   %imagesc(VelocityInBlock)
   %colorbar
   figure(6)
   surf(Z)
%--------------------------------------------------------------------------




