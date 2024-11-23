function [xs, ys, totalError] = calculateDistance( sensor1, sensor2, sensor3, t1, t2, t3, velocity, convergence_criteria )

dt1 = t2 - t1;
dt2 = t3 - t1;

c = velocity;

x1 = sensor1(1);
y1 = sensor1(2);

x2 = sensor2(1);
y2 = sensor2(2);

x3 = sensor3(1);
y3 = sensor3(2);

D1 = sqrt( (x2-x1)^2 + (y2-y1)^2 );
D2 = sqrt( (x3-x1)^2 + (y3-y1)^2 );

theta1 = acos( (x2-x1)/D1 );
theta3 = acos( (x3-x1)/D2 );

converged = false;

%initial theta guess
theta_current = 0;
i = 1;

d1_1 = [];
d1_2 = [];
error = [];
theta = [];
i = 1;

while theta_current <= 2*pi
    
    d1_1(i) = ( D1^2 - dt1^2 * c^2 ) / ( 2 * ( dt1*c + D1*cos(theta_current - theta1) ) );

    d1_2(i) = ( D2^2 - dt2^2 * c^2 ) / ( 2 * ( dt2*c + D2*cos(theta3 - theta_current) ) );
    
    error(i) = abs(d1_2(i) - d1_1(i));
    
    theta(i) = theta_current;
    
    theta_current = theta_current + 0.000001;
    
    i = i+1;
end


indexOfMinError = find( error == min(error), 1, 'first');
bestTheta = theta(indexOfMinError);
best_d1_1 = d1_1(indexOfMinError);
best_d1_2 = d1_2(indexOfMinError);
avg_d = (best_d1_1+best_d1_2)/2;

x1 = sensor1(1);
y1 = sensor1(2);

xs = x1 + avg_d*cos(bestTheta)
ys = y1 + avg_d*sin(bestTheta)

xs_real = 0.0;
ys_real = 0.0;

totalError = sqrt((xs_real - xs)^2 + (ys_real - ys)^2)

end
