%% Apexes.m
% Function to calculate the apexes of a track.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function ApexIndizes = Apexes(R)

k = 2;
ApexIndizes = [];

while k < length(R)-1
    
    if R(k-1) > R(k) && R(k+1) > R(k) && R(k) < 400
        ApexIndizes = [ApexIndizes k];
    end
    
    k = k+1;
end

end