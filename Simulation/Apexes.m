function ApexIndizes = Apexes(R)

k = 2;
ApexIndizes = [];

while k < length(R)-1
    
    if R(k-1) > R(k) && R(k+1) > R(k) && R(k) < 40
        ApexIndizes = [ApexIndizes k];
    end
    
    k = k+1;
end

end