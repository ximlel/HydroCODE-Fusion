function z = index(x,y)
%center=0.3;  %center of the interface seperated two fluid.
%A=0.02;  %amplitude
d=(0.35+0.02*cos(x*2.0*pi/5.93333)-y)*10.0;
if d<=0.0
	z=1.0;
elseif d>=1.0
	z=0.0;
else
	z=exp(log(1e-15)*d^8);
end
