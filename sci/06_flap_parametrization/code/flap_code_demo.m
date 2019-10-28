
close all
SDD = shape_dynamizer();
param = SF.parameters.full_non_linear.full_export;
SDD.flap_angle = -20;
SDD.flap_knee =  0.1;
SDD.x_hinge = 0.5;
[tx, tz] = SD.generate_coordinates(200, param);
[s,n] = SDD.T_curvilinear_transformation(tx, tz);
x_range = linspace(0, 1, 20);
y_range = linspace(-0.2, 0.2, 11);
[X, Y] = meshgrid(x_range, y_range);
[S, N] = SDD.T_curvilinear_transformation(X, Y);

% Plot deformed stuff!
figure(2)
surf(S, N, zeros(size(S)));
hold on
plot(s, n)
view(2)
axis([0 , 1.5, -0.75 ,0.75])

legend('Deformed Grid' , 'Deformed Profile', 'Location', 'SouthEast')
title('Flap as an Orthogonal Curvilinear Transformation - Merci Lame (1859)')
xlabel('eta')

set(gcf, 'PaperType', 'A5')
orient landscape
print('-dpdf', 'deformed_airfoil_on_grid.pdf')

% Plot un-deformed stuff!
figure(3)
surf(X, Y, zeros(size(S)));
hold on
plot(tx, tz)
view(2)
axis([0 , 1.5, -0.75 ,0.75])

atan2(N(3, 10) , S(3, 10)-SDD.x_hinge) * 180 / pi()

legend('Original Grid' , 'Original Profile', 'Location', 'SouthEast')
title('Flap as an Orthogonal Curvilinear Transformation - Merci Lame (1859)')
xlabel('eta')

set(gcf, 'PaperType', 'A5')
orient landscape
print('-dpdf', 'original_airfoil_on_grid.pdf')

disp('Transformation Accuracy Measures:')
cos_original = acos(dot([X(2, end ) - X(1, end) , Y(2, end ) - Y(1, end)] , [X(1, end-1 ) - X(1, end) , Y(1, end-1) - Y(1, end)]) ...
    / (norm([X(2, end ) - X(1, end) , Y(2, end ) - Y(1, end)]) * norm([X(1, end-1 ) - X(1, end) , Y(1, end-1) - Y(1, end)]))) * 180 / pi();

cos_deformed = acos(dot([S(2, end ) - S(1, end) , N(2, end ) - N(1, end)] , [S(1, end-1 ) - S(1, end) , N(1, end-1) - N(1, end)]) ...
    / (norm([S(2, end ) - S(1, end) , N(2, end ) - N(1, end)]) * norm([S(1, end-1 ) - S(1, end) , N(1, end-1) - N(1, end)]))) * 180 / pi();
disp(['    Shear between grid line versors     : ' num2str(cos_deformed - cos_original) ' (degrees, @ NE corner)']);

mesh_relative_compression = norm([S(1, end-1 ) - S(1, end) , N(1, end-1) - N(1, end)]) / norm([X(1, end-1 ) - X(1, end) , Y(1, end-1) - Y(1, end)]);
disp(['    Chorwise Cell Length preservation   : ' num2str(mesh_relative_compression * 100) '     (%      , @ NE corner)'])


hingeline_angle_deformed = atan2(N(6,end) - N(6, end-1) , S(6,end) - S(6,end-1)) * 180 / pi();
hingeline_angle_original= atan2(Y(6,end) - Y(6, end-1) , X(6,end) - X(6,end-1)) * 180 / pi();
disp(['    Original Longitudinal Mesh Angle    : ' num2str(hingeline_angle_original)]);
disp(['    Deformed Longitudinal Mesh Angle    : ' num2str(hingeline_angle_deformed)]);

% Make a reference hingeline vector
tx = linspace(0,1,1000);
tz = ones(size(tx)) * SDD.z_hinge;

% Now transform the hingeline
[ts,tn] = SDD.T_curvilinear_transformation(tx, tz);
tz2 = SDD.hingeline_deformation(tx);

% And plot
figure(2)
plot(tx, tz, tx, tz2, '--', ts, tn);
grid on
legend('Deformed Grid' , 'Deformed Profile', 'Original Hingeline' , 'Deformed Hingeline', 'Deformed and Renormalized Hingeline' , 'Location', 'SouthEast')
axis([-0.1 1.1 -0.6 0.6])

original_hingeline_length = sum(sqrt(diff(tx).^2 + (diff(tz)).^2));
deformed_hingeline_length = sum(sqrt(diff(ts).^2 + (diff(tn)).^2));
disp(['    Original Hingeline Length           : ' num2str(original_hingeline_length)])
disp(['    Deformed Hingeline Length           : ' num2str(deformed_hingeline_length)])


effective_hingeline_angle = atan2(tn(end) - tz(1), ts(end)-SDD.x_hinge) * 180 / pi();
disp(['    Effective hingeline angle           : ' num2str(effective_hingeline_angle)])
