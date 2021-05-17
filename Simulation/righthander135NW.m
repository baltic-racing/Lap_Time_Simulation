
function h = righthander135NW(x,y,r)
alpha = 1.25*pi:-pi/1000:0.5*pi;
xCoordinate = r * cos(alpha) + x;
yCoordinate = r * sin(alpha) + y;
h = [xCoordinate' yCoordinate'];
end