
function h = righthander90N(x,y,r)
alpha = pi:-pi/1000:0.5*pi;
xCoordinate = r * cos(alpha) + x;
yCoordinate = r * sin(alpha) + y;
h = [xCoordinate' yCoordinate'];
end