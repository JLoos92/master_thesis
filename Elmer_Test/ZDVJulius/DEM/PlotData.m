clear all
close all

d = load('BEDplusBump.xyz')
scatter(d(:,1),d(:,2),4,d(:,3))