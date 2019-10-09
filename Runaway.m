% Runway code for rendezvous_tmin
n=6;
x=[3115.66034733000;-104.839566850000;2814.27072303000;-1260.31722449000;-844.933173910000;1469.00192739000];
y=[-6442.03946922000;-3616.21292276000;-6560.49198928000;5130.94001982000;1921.55249890000;-6987.65249164000];
z=[0.00291745000000000;5931.68961103000;0.00418982000000000;4822.90029814000;6829.82848400000;0.00732756000000000];
xdot=[0.418370500000000;1.10036539000000;0.425565970000000;-1.82157690000000;-4.27048433000000;0.460030470000000];
ydot=[0.216240790000000;6.39563389000000;0.198809550000000;4.72943816000000;5.74646706000000;0.102254240000000];
zdot=[7.45264753000000;3.92012750000000;7.46259588000000;-5.47818080000000;-2.12905271000000;7.45493993000000];

[sequence,min_t]=rendezvous_tmin(6,x,y,z,xdot,ydot,zdot,2) %1- Curtis method %other than 1- Prussing method
           