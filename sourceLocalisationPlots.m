%Source Localisation Visualisation Script

source_location = [0.0, 0.0];

sensor1 = [0, 0.0090];
sensor2 = [ 0, 0.0040];
sensor3 = [0, 0.0340];
%lines for plotting


S0_input_start_times = [1.19E-09, 1.26E-09, 0.007555437;...
    0.012633409, 0.012850012, 0.014912133; 0.012633409, 0.012850012, 0.016716947];

A0_input_start_times = [1.14E-09, 1.22E-09, 0.007582028;...
    0.012493779, 0.012718204, 0.014750885; 0.012633409, 0.012850012, 0.016716947];

Comb_input_start_times = [1.19E-09, 1.26E-09, 0.007555437;...
    0.012633409, 0.012850012, 0.014912133; 0.012633409, 0.012850012, 0.016716947];


S0_input_max_times = [-2.82E-09, -2.52E-09, 0.00158193;...
    0.00038423, 0.00128702, 0.006780666; 0.00038423, 0.00128702, 0.006962753];

A0_input_max_times = [-1.25E-08, -1.16E-08, -6.18E-09;...
    -0.029126592, -0.026570633, -0.009908299; 0.029126592, 0.026570633, 0.009908299];

Comb_input_max_times = [-2.82E-09, -2.52E-09, 0.00158193;...
    3.84E-04, 0.00128702, 0.006780666; 3.84E-04, 0.00128702, 0.006962753];


sensor1 = [0, 0.0090];
sensor2 = [ 0, 0.0040];
sensor3 = [0, 0.0340];

figure
subplot(1,2,1)
scatter(source_location(1), source_location(2), 20, 'kd','filled')
hold all
scatter(sensor1(1), sensor1(2),15, 'k')
hold on


scatter(Comb_input_start_times(1,1), Comb_input_start_times(2,1),30, 'xr')
scatter(Comb_input_start_times(2,2), Comb_input_start_times(2,2),30, 'xg')
scatter(Comb_input_start_times(1,3), Comb_input_start_times(2,3),30, 'xb')

scatter(Comb_input_max_times(1,1), Comb_input_max_times(2,1),30, '+r')
scatter(Comb_input_max_times(2,2), Comb_input_max_times(2,2),30, '+g')
scatter(Comb_input_max_times(1,3), Comb_input_max_times(2,3),30, '+b')

legend("Signal Source", "Sensor", "T_0 Prediction, Extensional",...
    "T_0 Prediction, Longitudinal", "T_0 Prediction, Transverse", "T_{Max} Prediction, Extensional",...
    "T_{Max} Prediction, Longitudinal", "T_{Max} Prediction, Transverse",'AutoUpdate','off', 'location', 'southwest');

scatter(sensor2(1), sensor2(2),15, 'k')
scatter(sensor3(1), sensor3(2),15, 'k')
rectangle('Position',[-0.25 -0.3 0.5 0.6])
axis([-0.3 0.3 -0.35 0.35])
xlabel("x (m)")
ylabel("y (m)")

subplot(1,2,2)
scatter(source_location(1), source_location(2), 30, 'kd','filled')
hold all
scatter(sensor1(1), sensor1(2),25, 'k')
hold on


scatter(Comb_input_start_times(1,1), Comb_input_start_times(2,1),45, 'xr')
scatter(Comb_input_start_times(2,2), Comb_input_start_times(2,2),45, 'xg')
scatter(Comb_input_start_times(1,3), Comb_input_start_times(2,3),45, 'xb')

scatter(Comb_input_max_times(1,1), Comb_input_max_times(2,1),45, '+r')
scatter(Comb_input_max_times(2,2), Comb_input_max_times(2,2),45, '+g')
scatter(Comb_input_max_times(1,3), Comb_input_max_times(2,3),45, '+b')

legend("Signal Source", "Sensor", "T_0 Prediction, Extensional",...
    "T_0 Prediction, Longitudinal", "T_0 Prediction, Transverse", "T_{Max} Prediction, Extensional",...
    "T_{Max} Prediction, Longitudinal", "T_{Max} Prediction, Transverse",'AutoUpdate','off', 'location', 'southeast');

scatter(sensor2(1), sensor2(2),25, 'k')
scatter(sensor3(1), sensor3(2),25, 'k')
rectangle('Position',[-0.25 -0.3 0.5 0.6])
axis([-0.03 0.03 -0.040 0.040])
xlabel("x (m)")
ylabel("y (m)")




errorS0StartTimes = S0_input_start_times(3,:);
errorS0MaxAmp = S0_input_max_times(3,:);
errorsS0 = [errorS0MaxAmp; errorS0StartTimes];

errorA0StartTimes = A0_input_start_times(3,:);
errorA0MaxAmp = A0_input_max_times(3,:);
errorsA0 = [errorA0MaxAmp; errorA0StartTimes];


errorCombStartTimes = Comb_input_start_times(3,:);
errorCombMaxAmp = Comb_input_max_times(3,:);
errorsComb = [errorCombMaxAmp; errorCombStartTimes];


figure

X = categorical({'C_E','C_L','C_T'});
X = reordercats(X,{'C_E','C_L','C_T'});

bar(X, errorsS0*1000)
ylabel("Error (mm)")
legend('T_{Max}','T_0', 'Location', 'northwest' )

figure

X = categorical({'C_E','C_L','C_T'});
X = reordercats(X,{'C_E','C_L','C_T'});

bar(X, errorsA0*1000)
ylabel("Error (mm)")
legend('T_{Max}','T_0', 'Location', 'northwest' )

figure

X = categorical({'C_E','C_L','C_T'});
X = reordercats(X,{'C_E','C_L','C_T'});

bar(X, errorsComb*1000)
ylabel("Error (mm)")
legend('T_{Max}','T_0', 'Location', 'northwest' )



