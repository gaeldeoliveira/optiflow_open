function AIRFOIL=calculate_structural_properties_airfoils_box(COORD,POINTS_OUTPUT,box)
% interpolate for 1000 points
COORD2=interpolate_1000_points(COORD,5000);

% normalize and rotate
COORD3=normalize_and_rotate(COORD2);

% calculate area (?? last argument ??)
AIRFOIL.AREA=polyarea([COORD3(:,1);COORD3(1,1)],[COORD3(:,2);COORD3(1,2)]);

%%%%%%%%% calculate thickness and camber
[THICKNESS_MAX,THICKNESS_MAX_POSITION,THICKNESS_LINE,CAMBER_MAX,CAMBER_MAX_POSITION,CAMBER_LINE]=calculate_camber_thickness(COORD3,100);
AIRFOIL.THICKNESS_MAX=THICKNESS_MAX;
AIRFOIL.THICKNESS_MAX_POSITION=THICKNESS_MAX_POSITION;
AIRFOIL.THICKNESS_LINE=THICKNESS_LINE;
AIRFOIL.CAMBER_MAX=CAMBER_MAX;
AIRFOIL.CAMBER_MAX_POSITION=CAMBER_MAX_POSITION;
AIRFOIL.CAMBER_LINE=CAMBER_LINE;

%%%% calculate bending properties
[XBEND,YBEND]=calculate_bending_parameters(COORD3,box);
AIRFOIL.XBEND=XBEND;
AIRFOIL.YBEND=YBEND;

%%%% calculate torsional stiffness 
DX=COORD3(2:end,1)-COORD3(1:end-1,1);
DY=COORD3(2:end,2)-COORD3(1:end-1,2);
SL=sum(sqrt(DX.^2+DY.^2));
AIRFOIL.Jskin_over_t=4*AIRFOIL.AREA.^2/SL;

%%%% calculate LE radius and TE angle
[R_LE,ANGLE_TE]=calculate_LEradius_TEangle(COORD);
AIRFOIL.R_LE=R_LE;
AIRFOIL.ANGLE_TE=ANGLE_TE;

COORD4=interpolate_1000_points(COORD3,POINTS_OUTPUT);
AIRFOIL.COORD=COORD4;

% % % % % figure
% % % % % hold on 
% % % % % plot(COORD3(:,1),COORD3(:,2),'b-')
% % % % % plot(AIRFOIL.THICKNESS_LINE(:,1),AIRFOIL.THICKNESS_LINE(:,2),'g-')
% % % % % plot(AIRFOIL.CAMBER_LINE(:,1),AIRFOIL.CAMBER_LINE(:,2),'r-')
% % % % % 
% % % % % daspect([1 1 1])
% % % % % grid on
end


function COORD2=interpolate_1000_points(COORD,N)
% this function converts airfoil coordinates to N evenly spaced coordinates
DX=COORD(2:end,1)-COORD(1:end-1,1);
DY=COORD(2:end,2)-COORD(1:end-1,2);
DS=sqrt(DX.^2+DY.^2);
S=[0;cumsum(DS)];
L=0:max(S)/N:max(S);

X=interp1(S,COORD(:,1),L,'cubic');
Y=interp1(S,COORD(:,2),L,'cubic');

COORD2(:,1)=X;
COORD2(:,2)=Y;
end


function COORD2=normalize_and_rotate(COORD)
% finding trailing edge
TE=(COORD(1,:)+COORD(end,:))/2;

% find distance of points to trailing edge
R=sqrt((COORD(:,1)-TE(1)).^2+(COORD(:,2)-TE(2)).^2);

% find leading edge
[RMAX,IRMAX]=max(R);
IRMAX=IRMAX(1);
RMAX=RMAX(1);
LE=COORD(IRMAX,:);

% LE at 0,0 and normalize (RMAX is chord)
COORDTEMP=[COORD(:,1)-LE(1), COORD(:,2)-LE(2)]/RMAX;

% calculate angle and rotate
ANG=atan2(TE(2)-LE(2),TE(1)-LE(1));
COORD2(:,1)=  COORDTEMP(:,1)*cos(ANG)+COORDTEMP(:,2)*sin(ANG);
COORD2(:,2)= -COORDTEMP(:,1)*sin(ANG)+COORDTEMP(:,2)*cos(ANG);
end


function [THICKNESS_MAX,THICKNESS_MAX_POSITION,THICKNESS_LINE,CAMBER_MAX,CAMBER_MAX_POSITION,CAMBER_LINE]=calculate_camber_thickness(COORD,M)
% assumes normalized and rotated coordinates LE=0,0 and TE=1,0

% find leading edge
P=sum(COORD.^2,2);
[TEMP,ILE]=min(P);
ILE=ILE(1);

XI=[0:1/M:1]';

% UPPERSURFACE
YUPPER=interp1(COORD(1:ILE,1),COORD(1:ILE,2),XI,'cubic');
% LOWERSURFACE
YLOWER=interp1(COORD(ILE:end,1),COORD(ILE:end,2),XI,'cubic');

THICKNESS_LINE=[XI,YUPPER-YLOWER];
[THICKNESS_MAX,IMAX]=max(THICKNESS_LINE(:,2));
THICKNESS_MAX=THICKNESS_MAX(1);
IMAX=IMAX(1);
THICKNESS_MAX_POSITION=THICKNESS_LINE(IMAX(1),1);

CAMBER_LINE=[XI,(YUPPER+YLOWER)/2];
[CAMBER_MAX,IMAX]=max(CAMBER_LINE(:,2));
CAMBER_MAX=CAMBER_MAX(1);
IMAX=IMAX(1);
CAMBER_MAX_POSITION=CAMBER_LINE(IMAX(1),1);
end


function [R_LE,ANGLE_TE]=calculate_LEradius_TEangle(COORD)
%%% assumes normalized and rotated coordinates LE=0,0 and TE=1,0

%%% find leading edge
P=sum(COORD.^2,2);
[TEMP,ILE]=min(P);
ILE=ILE(1);

%%%% surface coord
DX=COORD(2:end,1)-COORD(1:end-1,1);
DY=COORD(2:end,2)-COORD(1:end-1,2);
DS=sqrt(DX.^2+DY.^2);
S=[0;cumsum(DS)];

%%%%%%%%%%%%%%% calculate curvature
X=COORD(:,1);
Y=COORD(:,2);

XL=gradient(X)./gradient(S);
XLL=gradient(XL)./gradient(S);

YL=gradient(Y)./gradient(S);
YLL=gradient(YL)./gradient(S);

KCURV=abs(XL.*YLL-YL.*XLL)./(XL.^2+YL.^2).^(3/2);

R_LE=1/KCURV(ILE);

%%%%%%%%%% calcualte TE angle
A1=atan(-(Y(2)-Y(1))/(X(2)-X(1)));
A2=atan((Y(end-1)-Y(end))/(X(end-1)-X(end)));

ANGLE_TE=(A1+A2)*180/pi;
end

function [XBEND,YBEND]=calculate_bending_parameters(COORD,box)
%%% assumes normalized and rotated coordinates LE=0,0 and TE=1,0

%%% temp number of slices
M=5000;
%%% find leading edge
P=sum(COORD.^2,2);
[TEMP,ILE]=min(P);
ILE=ILE(1);

XI=[0:1/M:1]';
DX=XI(2)-XI(1);

%% box definition, front and back of box location on chord
box1 = find(XI == box.front);
box2 = find(XI == box.back);

XI = XI(box1:box2);

%%% UPPERSURFACE
YUPPER=interp1(COORD(1:ILE,1),COORD(1:ILE,2),XI,'cubic');
%%% LOWERSURFACE
YLOWER=interp1(COORD(ILE:end,1),COORD(ILE:end,2),XI,'cubic');

DY=YUPPER-YLOWER;

% find midpoint of slices
XSLICES=(XI(2:end)+XI(1:end-1))/2;  
YSLICES=(YUPPER(2:end)+YUPPER(1:end-1)+YLOWER(2:end)+YLOWER(1:end-1))/4;

AVERAGE_HEIGHT_SLICES=(DY(2:end)+DY(1:end-1))/2;

AREA_SLICES=DX*AVERAGE_HEIGHT_SLICES;



%% calculate Xc and Yc centroid
XBEND.XC=sum(XSLICES.*AREA_SLICES)./sum(AREA_SLICES);
YBEND.YC=sum(YSLICES.*AREA_SLICES)./sum(AREA_SLICES);

%% maximum and minimum distances to centroids
XBEND.MAX_XminusXc=1-XBEND.XC;
XBEND.MIN_XminusXc=-XBEND.XC;
YBEND.MAX_YminusYc=mean(max(COORD(:,2)-YBEND.YC));
YBEND.MIN_YminusYc=mean(min(COORD(:,2)-YBEND.YC));

%% Define COORD into box coordinates
half = round(length(COORD)/2);      % find half of points
COORD_up = COORD(1:half,:);         % split up two sides of airfoil to find both 
COORD_low = COORD(half+1:end,:);

[~,i_front_up] = min(abs(COORD_up(:,1)-box.front));         % find positions of box webs
[~,i_front_low] = min(abs(COORD_low(:,1)-box.front));
[~,i_back_up] = min(abs(COORD_up(:,1)-box.back));
[~,i_back_low] = min(abs(COORD_low(:,1)-box.back));

%%%% calculate values of Iyy
XBEND.SOLID_Iyy=sum(AVERAGE_HEIGHT_SLICES.*(-(XI(1:end-1)-XBEND.XC).^3+(XI(2:end)-XBEND.XC).^3))/3; % only solid box
    DXL=COORD(2:end,1)-COORD(1:end-1,1);
    DYL=COORD(2:end,2)-COORD(1:end-1,2);
    DS=sqrt(DXL.^2+DYL.^2);
Iyy_over_tsurf_element = DS.*(((COORD(2:end,1)-XBEND.XC).^2+(COORD(1:end-1,1)-XBEND.XC).^2)/2); % find Iyy/t for each element

XBEND.Iyy_over_tsurf=sum(Iyy_over_tsurf_element);                               % sum to get total Iyy/t of airfoil
% Iyy/t parameters of box
XBEND.Iyy_box_skin = sum(Iyy_over_tsurf_element(i_back_up:i_front_up))...   
    + sum(Iyy_over_tsurf_element(half+1+i_front_low:half+1+i_back_low));
XBEND.Iyy_box_webs = box.web_t*DY(1)*(box.front-XBEND.XC)^2 + box.web_t*DY(end)*(box.back-XBEND.XC)^2;  % term with 1/12 is left out (thinwalled)

%%%% calculate values of Ixx
YBEND.SOLID_Ixx=  sum(DX.*(((YUPPER(1:end-1)-YBEND.YC).^3+(YUPPER(2:end)-YBEND.YC).^3)/2 ...
    -((YLOWER(1:end-1)-YBEND.YC).^3+(YLOWER(2:end)-YBEND.YC).^3)/2)/3);                     % only solid box

    DX=COORD(2:end,1)-COORD(1:end-1,1);
    DY=COORD(2:end,2)-COORD(1:end-1,2);
    DS=sqrt(DX.^2+DY.^2);
Ixx_over_tsurf_element = DS.*(((COORD(2:end,2)-YBEND.YC).^2+(COORD(1:end-1,2)-YBEND.YC).^2)/2);

YBEND.Ixx_over_tsurf = sum(Ixx_over_tsurf_element);
YBEND.Ixx_box_skin = sum(Ixx_over_tsurf_element(i_back_up:i_front_up))...
    + sum(Ixx_over_tsurf_element(half+1+i_front_low:half+1+i_back_low));
YBEND.Ixx_box_webs = 1/12*box.web_t*DY(1)^3 + box.web_t*DY(1)*((YUPPER(1)+YLOWER(1))/2-YBEND.YC)^2 ...
    + 1/12*box.web_t*DY(end)^3 + box.web_t*DY(end)*((YUPPER(end)+YLOWER(end))/2-YBEND.YC)^2;

end











