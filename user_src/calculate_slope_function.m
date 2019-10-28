function [SCMAX,ANGLE_SCMAX,SC_ANGLE_IN]=calculate_slope_function(ap,RANGE_ANGLE_OPERATION,ANGLE_IN)
% inputs:   ap                      = simulation results
%           RANGE_ANGLE_OPERATION   = range of angle of attack operation of VAWT (depending on TSR)
%           ANGLE_IN                = angle at which SC has to be specified (used for coupling between 2 sims)
%
% outputs:  SCMAX                   = maximum lift slope over drag
%           ANGLE_SCMAX             = angle at which SCMAX occurs
%           SC_ANGLE_IN             = maximum lift slope at specified angle (used for coupling between 2 sims)          

%%%%% weighted input values
DIST=[0.392208035249415   0.088932633017998   0.060492358783630   0.051692198661402   0.044899149911516 ...
    0.047310683677405   0.049033417456065   0.053759198899472   0.066593195839614   0.145079128503482];


AMIN=min(ap.raw_data.alpha);
AMAX=max(ap.raw_data.alpha);
%[MAXCL,IMAXCL]=max(ap.raw_data.cl);
%[MINCL,IMINCL]=min(ap.raw_data.cl);

%% interpolate data into a regular angle distribution
DANGL=0.01;
ANG=(AMIN):DANGL:(AMAX);
CLT=interp1(ap.raw_data.alpha,ap.raw_data.cl,ANG,'linear');
CDT=interp1(ap.raw_data.alpha,ap.raw_data.cd,ANG,'linear');

%% determine DCL_DALPHA
DCDALPHAT=gradient(CLT)./gradient(ANG);

%% add fake data to control filters
DCDALPHA=[-5000 DCDALPHAT  -5000];
CD=[50000 CDT 50000];
ANG=[AMIN-DANGL ANG AMAX+DANGL];

NELEMENTS=round(RANGE_ANGLE_OPERATION/DANGL);
FILTER_DCL=ones(1,NELEMENTS)/NELEMENTS;

%%% weighted values
IND1=(1:1:NELEMENTS)/NELEMENTS*10;
IND1=floor(IND1)+1;
IND1(IND1>10)=10;
FILTER_CD=DIST(IND1);
FILTER_CD=FILTER_CD/sum(FILTER_CD);


DCL_DALPHA_AVERAGE_RANGE=filter2(FILTER_DCL,DCDALPHA);
DCD_AVERAGE_RANGE=filter2(FILTER_CD,CD);
DCLCD_AVERAGE=DCL_DALPHA_AVERAGE_RANGE./(DCD_AVERAGE_RANGE);
[DCLCD_AVERAGEMAX,I1]=max(DCLCD_AVERAGE);
SCMAX=DCLCD_AVERAGEMAX(1);
ANGLE_SCMAX=ANG(I1(1));
SC_ANGLE_IN=interp1(ANG,DCLCD_AVERAGE,ANGLE_IN);

end