
function h = lefthander90(x,y,r)
alpha = -0.5*pi:pi/1000:0;
xCoordinate = r * cos(alpha) + x;
yCoordinate = r * sin(alpha) + y;
h = [xCoordinate' yCoordinate'];
end