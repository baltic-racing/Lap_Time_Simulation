function R = CurvatureRadius(X,Y)
K = 0;

for k = 2:length(X)-1
%     K(k) = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
%         sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));
    
    K(k) = 2*abs((X(k)-X(k-1)).*(Y(k+1)-Y(k-1))-(X(k+1)-X(k-1)).*(Y(k)-Y(k-1))) ./ ...
        sqrt(((X(k)-X(k-1)).^2+(Y(k)-Y(k-1)).^2)*((X(k+1)-X(k-1)).^2+(Y(k+1)-Y(k-1)).^2)*((X(k+1)-X(k)).^2+(Y(k+1)-Y(k)).^2));
end

R = 1./K;   % [m] Radius an jeder Stelle

end