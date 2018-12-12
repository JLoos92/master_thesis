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
xv=(-LX/2:dx:LX/2)+GLx;
yv=(-LY/2:dy:LY/2);
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
display(str)
str=strcat('Coefficients for V2 as function of height: ',num2str(coeffV2));
display(str)

figure(33)
subplot(1,3,1)
plot(x,z,'k.');hold on
plot(x(backboundary1),z(backboundary1),'gx')
plot(x(backboundary2),z(backboundary2),'rx')
grid on
subplot(1,3,2)
plot(V1(backboundary1),Height(backboundary1),'g-x');hold on;
plot(V1(backboundary2),Height(backboundary2),'r--')
plot(meanV1back,meanHeightback,'b--')
plot(polyval(coeffV1,meanHeightback),meanHeightback,'ko')
grid on
subplot(1,3,3)
plot(V2(backboundary1),Height(backboundary1),'g-x');hold on;
plot(V2(backboundary2),Height(backboundary2),'r--')
plot(meanV2back,meanHeightback,'b--')
plot(polyval(coeffV2,meanHeightback),meanHeightback,'ko')
grid on

%Extrusion

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
save('ZB.xyz','ZBouts','-ASCII')
save('ZS.xyz','ZSouts','-ASCII')
save('BED.xyz','BEDouts','-ASCII')

%Write 4 point Mesh.geo file for GMSS
%It will include 4 boundaries.
fileID = fopen('Mesh.geo','w');
corners(1,:) = [max(yv)-dist2boundary,min(xv)+dist2boundary]; %rightback looking downstream
corners(2,:) = [min(yv)+dist2boundary,min(xv)+dist2boundary];  %leftback
corners(3,:) = [min(yv)+dist2boundary,max(xv)-dist2boundary]; %leftfront
corners(4,:) = [max(yv)-dist2boundary,max(xv)-dist2boundary]; %rightfront

fileID = fopen('Mesh.geo','w');
fprintf(fileID,'Mesh.Algorithm=5;\n');
for k=1:4
   str = strcat('Point(',num2str(k),')={',num2str(corners(k,2)),',',num2str(corners(k,1)),',0.0,',num2str(meshresolution),'};\n');
   fprintf(fileID,str);
end

fprintf(fileID,'Line(1)={1,2};\n');
fprintf(fileID,'Line(2)={2,3};\n');
fprintf(fileID,'Line(3)={3,4};\n');
fprintf(fileID,'Line(4)={4,1};\n');
fprintf(fileID,'Line Loop(5)={1,2,3,4};\n')
fprintf(fileID,'Physical Line(6)={1};\n')
fprintf(fileID,'Physical Line(7)={2};\n')
fprintf(fileID,'Physical Line(8)={3};\n')
fprintf(fileID,'Physical Line(9)={4};\n')
fprintf(fileID,'Plane Surface(10)={5};\n')
fprintf(fileID,'Physical Surface(11)={10};\n')
fprintf(fileID,'Transfinite Surface {10};\n')
fprintf(fileID,'Recombine Surface {10};')
fclose(fileID)

str=strcat('$yzero=',num2str(min(yv)));
display(str)
str=strcat('$xzero=',num2str(min(xv)));
display(str)
[w h]=size(BED);
str=strcat('$llx=',num2str((w-1)*dx));
display(str)
str=strcat('$lly=',num2str((h-1)*dy));
display(str)
str=strcat('$nny=',num2str(h));
display(str)
str=strcat('$nnx=',num2str(w));
display(str)
str='$nanvalue = -9999.0';
display(str)

%Indices for areas of interest.
[val valindback] = min(abs(x-min(xv)));
backboundary = find(x==x(valindback));
[val valindfront] = min(abs(x-max(xv)));
frontboundary = find(x==x(valindfront));
indbox = find(x>x(valindback) & x<x(valindfront));

%figure(1)
%subplot(1,3,1)
%plot(x(indnonzeroBed)/1000,Bed,'k--')
%hold on
%plot(xv/1000,Bedi,'kx')
%plot(x(indnonzeroZs)/1000,Zs,'b--')
%plot(xv/1000,Zsi,'bx')
%plot(x(indnonzeroZb)/1000,Zb,'y--')
%plot(xv/1000,Zbi,'yx')
%xlabel('distance (km)')
%ylabel('height (m)')
%subplot(1,3,2)
%scatter(x/1000,z,5,V1)
%colorbar;
%subplot(1,3,3)
%scatter(x/1000,z,5,V2)
%colorbar

%figure(2)
%subplot(2,1,1)
%plot(x,z,'k.');hold on
%plot(x(indbox),z(indbox),'g.')

%V1i=griddata(x(indbox),y(indbox),V1(indbox),xv,y(indbox)

%figure(3)
%imagesc(yv/1000,xv/1000,ZS);set(gca(),'YDir','normal')
%hold on
%plot(corners(1,1)/1000,corners(1,2)/1000,'kx')
%plot(corners(2,1)/1000,corners(2,2)/1000,'rx')
%plot(corners(3,1)/1000,corners(3,2)/1000,'yx')
%plot(corners(4,1)/1000,corners(4,2)/1000,'mx')
%figure(4)
% Bump function with sigma and amplitude

M = zeros(251,51)
N = 1; % number of bumps
sigma = 1000;% std (width) of Gauss 
maxAmplitude = 200; % maximum height

  % random location of bumps or xc, yc for grounding line position
    xc = 0
    yc = GLx


sigma_x = 500;
sigma_y = 1000;

theta = 2*pi
    a = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
    b = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
    c = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);

    Z_new = maxAmplitude*exp( - (a*(X-xc).^2 + 2*b*(X-xc).*(Y-yc) + c*(Y-yc).^2));
 
    % Gauss function
    exponent = ((X-xc).^2 + (Y-yc).^2)./(sigma^2);
    amplitude = maxAmplitude; 
    
    % add Gauss to the matrix M
    M = M + amplitude*exp(-exponent);


% add Gauss to topography
Zero = zeros(251,51)
BED_B = BED + Z_new;
ZB_B = ZB + [Z_new(1:125,:);Zero(126:251,:)];
%ZB_B = [ZB_B, ZB_new(125:251,:)]

%figure(5)
%surf(BED_B)
%shading interp

% safe new bedrock topo for elmer
BEDout_1 =[Y(:), X(:), BED_B(:)];
ZBout_1 =[Y(:), X(:), ZB_B(:)];

BED_bump600 = sortrows(BEDout_1,1);
ZB_bump600 = sortrows(ZBout_1,1);

%save('/DEM/BEDbump200_5001000_0.xyz','BED_bump600','-ASCII');
%save('/DEM/ZBbump200_5001000_0.xyz','ZB_bump600','-ASCII');
NEW = load ('BED_bump300200020000_20000.xyz')

% figure for all

figure(6)
surf (ZB_B)
hold on
%surf (ZB)
surf(BED_B)
surf (ZS)

