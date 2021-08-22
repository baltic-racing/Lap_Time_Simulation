function ApexIndizes = Apexes(R)
% Function to calculate the apexes of a track

k = 2;
ApexIndizes = [];

while k < length(R)-1
    
    if R(k-1) > R(k) && R(k+1) > R(k) && R(k) < 40
        ApexIndizes = [ApexIndizes k];
    end
    
    k = k+1;
end

end