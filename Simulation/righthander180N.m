
function h = righthander180N(x,y,r)
alpha = pi:-pi/1000:0;
xCoordinate = r * cos(alpha) + x;
yCoordinate = r * sin(alpha) + y;
h = [xCoordinate' yCoordinate'];
end