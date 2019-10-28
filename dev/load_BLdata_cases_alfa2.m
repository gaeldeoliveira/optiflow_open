
%% TN 2133
% Make bl_run object to store data
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 48.768;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 2680914.886;       % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'TN2133';          % [string] - Identifier for data origin
data_table             = [0.304800000000000,32.4224687200000,0.000812800000000000,1.35000000000000;0.457200000000000,34.0678601900000,0.000939800000000000,1.38000000000000;0.609600000000000,35.6373640700000,0.00106680000000000,1.38000000000000;0.762000000000000,37.0123096600000,0.00119380000000000,1.34000000000000;0.914400000000000,38.4618471500000,0.00121920000000000,1.40000000000000;1.06680000000000,39.8587042700000,0.00139700000000000,1.38000000000000;1.37160000000000,42.4029276700000,0.00137160000000000,1.39000000000000;1.67640000000000,45.2518689800000,0.00162560000000000,1.36000000000000;1.98120000000000,46.8528265100000,0.00182880000000000,1.35000000000000;2.28600000000000,47.8324340500000,0.00236220000000000,1.37000000000000;2.59080000000000,47.7826863100000,0.00251460000000000,1.37000000000000;2.89560000000000,47.7328867300000,0.00312420000000000,1.38000000000000;3.20040000000000,48.1792299300000,0.00337820000000000,1.35000000000000;3.50520000000000,48.0804006300000,0.00414020000000000,1.36000000000000;3.81000000000000,48.2039055900000,0.00411480000000000,1.36000000000000;4.11480000000000,48.2778569100000,0.00426720000000000,1.39000000000000;4.41960000000000,48.5725362900000,0.00487680000000000,1.33000000000000;4.72440000000000,48.4745088700000,0.00528320000000000,1.38000000000000;5.02920000000000,48.5480479900000,0.00574040000000000,1.34000000000000;5.33400000000000,48.7680000000000,0.00571500000000000,1.35000000000000;5.48640000000000,48.4745088700000,0.00581660000000000,1.37000000000000;5.63880000000000,47.9317746200000,0.00662940000000000,1.31000000000000;5.79120000000000,46.9542396700000,0.00716280000000000,1.37000000000000;5.94360000000000,46.0076391800000,0.00779780000000000,1.33000000000000;6.09600000000000,45.0147396500000,0.00810260000000000,1.40000000000000;6.24840000000000,43.9724048800000,0.00906780000000000,1.45000000000000;6.40080000000000,42.9878232700000,0.00990600000000000,1.49000000000000;6.55320000000000,41.9518198600000,0.0112522000000000,1.47000000000000;6.70560000000000,40.7147089300000,0.0127254000000000,1.54000000000000;6.85800000000000,39.5892844800000,0.0157480000000000,1.60000000000000;7.01040000000000,38.5544892300000,0.0167640000000000,1.65000000000000;7.16280000000000,37.4276528600000,0.0180340000000000,1.75000000000000;7.31520000000000,36.4294022200000,0.0218440000000000,1.87000000000000;7.46760000000000,35.4701300900000,0.0241300000000000,1.99000000000000;7.62000000000000,34.7247337900000,0.0294640000000000,2.22000000000000;7.74192000000000,34.3390445000000,0.0302260000000000,2.39000000000000;7.85469600000000,34.0727466300000,0.0345440000000000,2.80000000000000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);

%% Bradhaw and Ferriss (Run A)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 36.576;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 2482991.333;       % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'RM3575A';         % [string] - Identifier for data origin
data_table             = [0.609600000000000,38.2914212100000,0.00140208000000000,1.46000000000000;0.762000000000000,36.7220125600000,0.00189738000000000,1.41100000000000;0.914400000000000,35.0061174500000,0.00264668000000000,1.45100000000000;1.06680000000000,33.7610547400000,0.00322834000000000,1.50000000000000;1.21920000000000,32.6531715100000,0.00356870000000000,1.50000000000000;1.52400000000000,30.9061789300000,0.00473710000000000,1.46500000000000;1.82880000000000,29.5111886900000,0.00611886000000000,1.47200000000000;2.13360000000000,28.3316477700000,0.00739648000000000,1.49500000000000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);

%% Bradhaw and Ferriss (Run B)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 36.576;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 2531280.771;       % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'RM3575B';         % [string] - Identifier for data origin
data_table             = [0.609600000000000,38.2914212100000,0.00224536000000000,1.40700000000000;0.762000000000000,36.7220125600000,0.00284480000000000,1.40000000000000;0.914400000000000,35.0061174500000,0.00347980000000000,1.46200000000000;1.06680000000000,33.7610547400000,0.00426720000000000,1.48800000000000;1.21920000000000,32.6531715100000,0.00495300000000000,1.46700000000000;1.52400000000000,30.9061789300000,0.00647700000000000,1.47100000000000;1.82880000000000,29.5111886900000,0.00800100000000000,1.48800000000000;2.13360000000000,28.3316477700000,0.00932180000000000,1.51500000000000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);

%% Bradhaw and Ferriss (Run C)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 36.576;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 2520303.951;       % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'RM3575C';         % [string] - Identifier for data origin
data_table             = [0.609600000000000,38.2914212100000,0.00330200000000000,1.36100000000000;0.762000000000000,36.7220125600000,0.00391160000000000,1.40300000000000;0.914400000000000,35.0061174500000,0.00505460000000000,1.37000000000000;1.06680000000000,33.7610547400000,0.00576580000000000,1.46600000000000;1.21920000000000,32.6531715100000,0.00657860000000000,1.45100000000000;1.52400000000000,30.9061789300000,0.00835660000000000,1.47700000000000;1.82880000000000,29.5111886900000,0.0100330000000000,1.53200000000000;2.13360000000000,28.3316477700000,0.0114554000000000,1.54300000000000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);

%% Aubertine and Eaton (2005, LDA)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 20.48;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 1330479.363;       % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'AUBERTINE2005';         % [string] - Identifier for data origin
data_table             = [-0.330000000000000,20.4800000000000,0.00251000000000000,1.34000000000000;0,20.6200000000000,0.00217000000000000,1.29000000000000;0.250000000000000,19.8700000000000,0.00310000000000000,1.39000000000000;0.330000000000000,19.6000000000000,0.00332000000000000,1.39000000000000;0.500000000000000,19.2000000000000,0.00378000000000000,1.42000000000000;0.670000000000000,18.8000000000000,0.00424000000000000,1.46000000000000;0.750000000000000,18.7100000000000,0.00460000000000000,1.48000000000000;1,18.2300000000000,0.00566000000000000,1.56000000000000;1.33000000000000,18.2600000000000,0.00537000000000000,1.47000000000000;1.67000000000000,17.5400000000000,0.00526000000000000,1.45000000000000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);

%% Skare and Krogstad (1994, LDA)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 22.35;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 1457543.471;      % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'SKARE1994';         % [string] - Identifier for data origin
data_table             = [3,22.3500000000000,0.0175600000000000,1.79300000000000;3.20000000000000,21.9200000000000,0.0201400000000000,1.84000000000000;3.40000000000000,21.2100000000000,0.0226000000000000,1.90100000000000;3.60000000000000,20.5300000000000,0.0248400000000000,1.93600000000000;3.80000000000000,20.1000000000000,0.0270300000000000,1.95700000000000;4,19.8100000000000,0.0301300000000000,2.00600000000000;4.20000000000000,19.4200000000000,0.0325400000000000,1.99900000000000;4.40000000000000,19.3800000000000,0.0348300000000000,1.98900000000000;4.60000000000000,18.8400000000000,0.0373500000000000,1.99800000000000;4.80000000000000,18.6700000000000,0.0389800000000000,1.99400000000000;5,18.3000000000000,0.0429800000000000,1.99800000000000;5.20000000000000,18.0400000000000,0.0457800000000000,1.98600000000000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);

%% Kays and Moffat A (Stanford, 1972)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 9.512808;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 6.44E+05;      % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'ANDERSEN1972A';         % [string] - Identifier for data origin
data_table             = [0.0508000000000000,9.51280800000000,0.000833120000000000,1.53658536600000;0.254000000000000,9.48842400000000,0.00132334000000000,1.48176583500000;0.558800000000000,9.49147200000000,0.00195072000000000,1.47786458300000;0.863600000000000,9.46404000000000,0.00248158000000000,1.44524053200000;1.16840000000000,9.46099200000000,0.00300990000000000,1.44219409300000;1.47320000000000,9.48842400000000,0.00349250000000000,1.43127272700000;1.77800000000000,9.46099200000000,0.00402336000000000,1.42361111100000;2.08280000000000,9.47013600000000,0.00455676000000000,1.41192865100000;2.28600000000000,9.46708800000000,0.00483108000000000,1.40378548900000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);

%% Kays and Moffat B (Stanford, 1972)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 9.049512;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 585385.4356;      % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'ANDERSEN1972B';         % [string] - Identifier for data origin
data_table             = [0.0508000000000000,9.04951200000000,0.000906780000000000,1.55742296900000;0.254000000000000,8.20826400000000,0.00182372000000000,1.58495821700000;0.558800000000000,7.48893600000000,0.00314452000000000,1.56381260100000;0.863600000000000,7.08964800000000,0.00431800000000000,1.53882352900000;1.16840000000000,6.78789600000000,0.00552704000000000,1.54227941200000;1.47320000000000,6.51052800000000,0.00658368000000000,1.51774691400000;1.77800000000000,6.32155200000000,0.00774446000000000,1.51197113800000;2.08280000000000,6.20268000000000,0.00868426000000000,1.49400409500000;2.28600000000000,6.11124000000000,0.00928624000000000,1.49042669600000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);


%% Kays and Moffat C (Stanford, 1972)
bl_run          = experimental_bl_run();
% Fill it in
bl_run.U_ref           = 8.878824;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
bl_run.U_ref_over_nu   = 593607.8642;      % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
bl_run.source          = 'ANDERSEN1972C';         % [string] - Identifier for data origin
data_table             = [0.0508000000000000,8.87882400000000,0.000911860000000000,1.56267409500000;0.254000000000000,7.83336000000000,0.00195580000000000,1.60649350600000;0.558800000000000,6.80008800000000,0.00379984000000000,1.65641711200000;0.863600000000000,6.25754400000000,0.00557530000000000,1.65284738000000;1.16840000000000,5.91007200000000,0.00720344000000000,1.64421720700000;1.47320000000000,5.65404000000000,0.00865124000000000,1.61391661800000;1.77800000000000,5.47725600000000,0.0100355400000000,1.59225512500000;2.08280000000000,5.30961600000000,0.0116611400000000,1.58200827700000;2.28600000000000,5.21208000000000,0.0131419600000000,1.59335137200000];
bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% Make experimental case
EC = experimental_case_bl_run(EDB);
% Load data into case
EC.load_bl_run(bl_run)
% Add case to database
result_flag = add_experimental_case(EDB, EC);



%% New case template (uncomment once!)
% bl_run          = experimental_bl_run();
% % Fill it in
% bl_run.U_ref           = ;            % [m/s   ] - Reference tunnel wind speed (may have been adjusted accross runs to correct for nu changes with temperature)
% bl_run.U_ref_over_nu   = ;       % [m^(-1)] - Reference Reynolds number per unit lenght (m, not ft! this is not an adimensional quantity)
% bl_run.source          = '';          % [string] - Identifier for data origin
% data_table             = ;
% bl_run.x               = data_table(:,1);   % [m     ] - Array of x positions at which B.L. was measured
% bl_run.Ue              = data_table(:,2);   % [m/s   ] - Array of edge velocity positions at measurement points 
% bl_run.theta           = data_table(:,3);   % [m     ] - Array of BL momentum displacements at measurement points 
% bl_run.H               = data_table(:,4);   % [adim. ] - Array of shape factors at measurement points 
% % Make experimental case
% EC = experimental_case_bl_run(EDB);
% % Load data into case
% EC.load_bl_run(bl_run)
% % Add case to database
% result_flag = add_experimental_case(EDB, EC);