
function sc = cost_function_StifnessFinal( ~ , experiments_results)
% Calculated weighted Cl/Cd of a polar at AOA/Weight combinations from
% alpha_i and w_i vectors

% Input extraction
% parameters can be ignored

% Extract aerodynamic polar from standard experiment input form
experiment_result   = experiments_results{1};
% parameters = experiment_result.x;
coordinates = experiment_result.coordinates;


ap = experiment_result.ap;
RANGE_ANGLE_OPERATION=14; % degrees range of operation
if strcmp(class(ap) , 'aerodynamic_polar')
    % If polar is valid
    
    % find the appropriate thickness (of the NACA00XX airfoil) interpolated between the simulations 
    
    [MAXCL,IMAXCL]=max(ap.raw_data.cl);
    [MINCL,IMINCL]=min(ap.raw_data.cl);
%    sc=(MAXCL-MINCL)./(ap.raw_data.alpha(IMAXCL(1))-ap.raw_data.alpha(IMINCL(end)));
    DANGL=0.1;
    ANG=min(ap.raw_data.alpha):DANGL:max(ap.raw_data.alpha);
    CL=interp1(ap.raw_data.alpha,ap.raw_data.cl,ANG,'linear');
    CD=interp1(ap.raw_data.alpha,ap.raw_data.cd,ANG,'linear');    
    DCDALPHA=gradient(CL)./gradient(ANG);
    NELEMENTS=round(RANGE_ANGLE_OPERATION/DANGL);
    FILTER_DCL=ones(1,NELEMENTS)/NELEMENTS;
    DC_DALPHA_AVERAGE_RANGE=filter2(FILTER_DCL,DCDALPHA);
    [MAX_DCL_DALPHA,IMAX]=max(DC_DALPHA_AVERAGE_RANGE);
    
    IS1=round(NELEMENTS/2);
    IS2=NELEMENTS-IS1;
    if size(IMAX)>1
    
       LL=length(IMAX);
       LL=round((LL+0.2)/2);
       IMAX=IMAX(LL);
    end    
%    t11=max((IMAX-IS1+1),1)
%    t12=
    DCLCD_AVERAGE=MAX_DCL_DALPHA/mean(CD(max(IMAX-IS1+1,1):min(IMAX+IS2-1,length(CD))));
    
    
    
    
else
    % If polar is invalid!
    DCLCD_AVERAGE = 0;
end











XCOR=interp1(1:1:length(coordinates.tx),coordinates.tx,1:.1:length(coordinates.tx));
YCOR=interp1(1:1:length(coordinates.tz),coordinates.tz,1:.1:length(coordinates.tz));
DL=sqrt((XCOR(2:end)-XCOR(1:end-1)).^2+((YCOR(2:end)-YCOR(1:end-1)).^2));

STIFNESS=sum(DL.*((YCOR(2:end)+YCOR(1:end-1)).^2)/4);

%%% calculate angle trailing edge

ANG1=atan((YCOR(2)-YCOR(1))/abs(XCOR(1)-XCOR(2)))-atan((YCOR(end-1)-YCOR(end))/abs(XCOR(end-1)-XCOR(end)));
ANG1=ANG1*180/pi;

MINANG=10;
CORRECTION= 1-(abs((min(0,ANG1-MINANG)))/MINANG)^2;

STIFNESS_CORRECTED_TE =STIFNESS*CORRECTION;


if DCLCD_AVERAGE < 10
 STIFNESS_CORRECTED_TE=0;
end

   % Display Information on Selection
%    disp(['     Tmax = ' num2str(round(MaxThickness*100)) '%    @ ' num2str(round(Location*100)) '% of the chord'])
   disp(['     STIFNESS = ' num2str(STIFNESS) ])
   disp(['     TE ANGLE = ' num2str(ANG1) ])
   disp(['     STIFNESS_CORRECTED_TE = ' num2str(STIFNESS_CORRECTED_TE) ])
 

   sc = STIFNESS_CORRECTED_TE;
end