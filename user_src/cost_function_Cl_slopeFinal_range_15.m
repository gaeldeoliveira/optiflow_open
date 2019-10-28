function sc = cost_function_Cl_slopeFinal_range_15(parameters , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors

%% Input extraction
% Extract parameters from parameters structure

% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{1};
ap = experiment_result.ap;
RANGE_ANGLE_OPERATION=15; % degrees range of operation
if strcmp(class(ap) , 'aerodynamic_polar')
    % If polar is valid
    
    % find the appropriate thickness (of the NACA00XX airfoil) interpolated between the simulations 
    
    [MAXCL,IMAXCL]=max(ap.raw_data.cl);
    [MINCL,IMINCL]=min(ap.raw_data.cl);
    sc=(MAXCL-MINCL)./(ap.raw_data.alpha(IMAXCL(1))-ap.raw_data.alpha(IMINCL(end)));
    DANGL=0.1;
    ANG=min(ap.raw_data.alpha):DANGL:max(ap.raw_data.alpha);
    CL=interp1(ap.raw_data.alpha,ap.raw_data.cl,ANG,'linear');
    CD=interp1(ap.raw_data.alpha,ap.raw_data.cl,ANG,'linear');    
    DCDALPHA=gradient(CL)./gradient(ANG);
    NELEMENTS=round(RANGE_ANGLE_OPERATION/DANGL);
    FILTER_DCL=ones(1,NELEMENTS)/NELEMENTS;
    DC_DALPHA_AVERAGE_RANGE=filter2(FILTER_DCL,DCDALPHA);
    [MAX_DCL_DALPHA,IMAX]=max(DC_DALPHA_AVERAGE_RANGE);
    
    sc=MAX_DCL_DALPHA;
    
else
    % If polar is invalid!
    sc = Inf;
end