
function h = righthander45W(x,y,r)
alpha = 1.5*pi:-pi/1000:1.25*pi;
xCoordinate = r * cos(alpha) + x;
yCoordinate = r * sin(alpha) + y;
h = [xCoordinate' yCoordinate'];
end