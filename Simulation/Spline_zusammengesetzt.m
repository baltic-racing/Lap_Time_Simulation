clear all
close all
clc

load('SplineData_1_70.mat')
TrackFIT = [Track];

y = load('SplineData_1_70_NEU.mat');
XX = y.Track(:,1);
YY = y.Track(:,2);

load('SplineData_70_100_NEU.mat')
TrackFIT = [TrackFIT; Track];

figure(1)
plot(TrackFIT(:,1),TrackFIT(:,2))
hold on
plot(XX,YY)
axis equal
