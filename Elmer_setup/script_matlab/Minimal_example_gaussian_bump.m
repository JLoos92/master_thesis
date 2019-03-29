%%%%%%% Minimal Example Bump %%%%%%%%%


clear all; close all
%MeshResolution
meshresolution=1000;
%Choose model domain smaller than data domain with this value.
dist2boundary=5000;
%Domain Width in m+5km to avoide interpolation issues at boundary; spacing;
LY = 50*1000;dy=1000;
%Domain Length in m; Center at GL
LX = 250*1000;dx=1000;
%Grounding line x in Mismip2D
GLx = 1054*1000;
%Vertical Layers for Velocity Interpolation
vertlayers=21;


%Load Input from Mismip2D
load('mismip_1a4_09.nh.csv')
Mismip2D = mismip_1a4_09_nh;
indnonzeroZs = find(Mismip2D(:,7)~=0);
indnonzeroZb = find(Mismip2D(:,9)~=0);
indnonzeroBed = find(Mismip2D(:,9)~=0);
Zs = Mismip2D(indnonzeroZs,7);
V1 = Mismip2D(:,20);
V2 = Mismip2D(:,21);
Depth = Mismip2D(:,11);
Height = Mismip2D(:,12);
Zb = Mismip2D(indnonzeroZb,9);
Bed = Mismip2D(indnonzeroBed,6);
x = Mismip2D(:,26);
z = Mismip2D(:,27);

%Get interpolated  2D DEM centered at GL and clipped to area of interest
xv=(-LX/2:dx:LX/2)+GLx;yv=(-LY/2:dy:LY/2);
Bedi = interp1(x(indnonzeroBed),Bed,xv,'PCHIP');
Zbi = interp1(x(indnonzeroZb),Zb,xv,'PCHIP');
Zsi = interp1(x(indnonzeroZs),Zs,xv,'PCHIP');

%Get interpolated velocities at back boundary
x2=x;
[val, valindback] = min(abs(x2-min(xv)));
backboundary1 = find(x2==x(valindback));
x2(backboundary1) = NaN;
[val2, valindback2] = min(abs(x2-min(xv)));
backboundary2 = find(x2==x(valindback2));

meanV1back = (val*V1(backboundary1)+val2*V1(backboundary2))/(val+val2);
meanV2back = (val*V2(backboundary1)+val2*V2(backboundary2))/(val+val2);
meanHeightback = (val*Height(backboundary1)+val2*Height(backboundary2))/(val+val2);


%Fit Polynomial, to be used at back boundary in Elmer
[coeffV1] = polyfit(meanHeightback,meanV1back,4);
[coeffV2] = polyfit(meanHeightback,meanV2back,2);
str=strcat('Coefficients for V1 as function of height: ',num2str(coeffV1));
display(str);
str=strcat('Coefficients for V2 as function of height: ',num2str(coeffV2));
display(str);

% Extrusion
 [X,Y] = meshgrid(yv,xv);
 
    for k=1:length(yv)
        ZB(:,k) = Zbi;
        BED(:,k)= Bedi;
        ZS(:,k) = Zsi;
    end
 
%Write DEM Files to be read into Elmer
ZBout =[Y(:), X(:), ZB(:)];
ZSout =[Y(:), X(:), ZS(:)];
BEDout =[Y(:), X(:), BED(:)];
ZBouts = sortrows(ZBout,1);
BEDouts = sortrows(BEDout,1);
ZSouts = sortrows(ZSout,1);
save('ZB.xyz','ZBouts','-ASCII');
save('ZS.xyz','ZSouts','-ASCII');
save('BED.xyz','BEDouts','-ASCII');

%Indices for areas of interest.
[val valindback] = min(abs(x-min(xv)));
backboundary = find(x==x(valindback));
[val valindfront] = min(abs(x-max(xv)));
frontboundary = find(x==x(valindfront));
indbox = find(x>x(valindback) & x<x(valindfront));


% Set bump at position x = GLx
new_t = normpdf(z(1000:4000),2000,1);

ind_glx = find(x==GLx)

sigma = 4000;% std (width) of Gauss 
maxAmplitude = 5000; % maximum height

 % Gauss function
    exponent = ((x-1100).^2 + (z-500).^2)./(2*sigma^2);
    amplitude = rand()*maxAmplitude;  
    % add Gauss to the matrix M
    M = x + amplitude*exp(-exponent);


% 2D Plot of iceshelf/sheet
figure(1)
subplot(2,1,1)
plot(M,z,'k.');hold on
%plot(x(indbox),z(indbox),'g.')
grid minor
