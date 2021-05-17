x(1) = 0;
y(1) = 0;
z(1) = 0;
s(1) = 0;
R(1) = 0;

for i = 2:300
    y(i) = 0;
    z(i) = 0;
    s(i) = s(i-1) + 0.25;
    x(i) = s(i);
    R(i) = 0;
end

Track.x_Track = x;
Track.y_Track = y;
Track.z_Track = z;
Track.s = s;
Track.R = R;

% save('AccelerationTrack', '-struct','Track');
save('AccelTrack.mat','Track')