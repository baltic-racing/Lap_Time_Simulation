
function h = righthander90S(x,y,r)
alpha = 2*pi:-pi/1000:1.5*pi;
xCoordinate = r * cos(alpha) + x;
yCoordinate = r * sin(alpha) + y;
h = [xCoordinate' yCoordinate'];
end