close all
keepdb


x = [0 0;1 0; 1 0; 1 0;0 1;0 0; NaN NaN ; NaN NaN ; -1 -1 ; NaN NaN ; NaN NaN ];
y = [1 1; 2 2; 4 3; 4 3];
P = polyline( x , y )


P.isclosed(1)

pplot(flip(P),'s-')

double(flip(P))


x = rand(1000,2)*8 - 2;
[p,y,d,c] = closestElement( P , x )
hplot3d( x , y , '.-r' )

%%

t = linspace(0,2*pi/3,11);
P = polyline( [ cos(t(:)) , sin(t(:)) ] );
P = P * [1 1] + [2 1];
P = double( P );
P = P + randn( size(P) )/20;
P = polyline( P );

d = 0.1;
PP = resample(P ,'e',d);
maxnorm( diff(PP.arclength{1}) - d )

pplot(P); axis equal
hplot( PP + [0 0] ,'.-r' )

%%

plotfun( @(d)NumberOfNodes( resample(P ,'e',d) ) , linspace(0,1,1000) );
maxnorm( diff(PP.arclength{1}) - d )

pplot(P); axis equal
hplot( PP + [0 0] ,'.-r' )





