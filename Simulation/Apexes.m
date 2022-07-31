%% Apexes.m
% Function to calculate the apexes of a track.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function ApexIndizes = Apexes(R)

% %% Old apex calculation algorithm
% k = 2;
% ApexIndizes = [];
% 
% while k < length(R)-1
%     
%     if R(k-1) > R(k) && R(k+1) > R(k) && R(k) < 300
%         ApexIndizes = [ApexIndizes k];
%     end
%     
%     k = k+1;
% end

%% New apex calculation algorithm
% ToDo: add parameters to track?
minPeakHeight = 0.0001; % 0.001
%MinPeakDistance = 100; % 200
MinPeakDistance = 30; % 200
MinPeakProminence = 0.005; % 0.001


% Ottobiano
minPeakHeight = 0.001; % 0.001
MinPeakDistance = 10; % 200 (Higher number less sensitive = less apexes, Minimum Distance from one peak (Apex) to another)
MinPeakProminence = 0.001; % 0.001

[~, ApexIndizes] = findpeaks(abs(1./R),'MinPeakDistance',MinPeakDistance,'MinPeakHeight',minPeakHeight,'MinPeakProminence',MinPeakProminence);

end