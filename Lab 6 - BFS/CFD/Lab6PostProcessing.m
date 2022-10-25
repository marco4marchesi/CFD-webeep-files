%{
this script runs only on my pc, to run it also on your pc you have to
download the repo at link:

https://github.com/marco4marchesi/progetto_CFD.git

and search for the matlab functions needed.


%}

clear; close all; clc;

addpath('C:\Users\marco\Desktop\UNI\2 MAGISTRALE\CFD\CFD PROJECT\progetto_CFD\matlabFunctions\')
addpath("plotOverLines\")

%% import data from simulations:

all_logs = folderLogExtractor("plotOverLines","raw");
len = length(all_logs.plotOverLine_x0.arc_length);
plotVec = ones(len,1);
ymin1 = 0;
ymax1 = 0.04;
ymin2 = -0.0267;
ymax2 = 0.04;
span1 = linspace(ymin1,ymax1,len);
span2 = linspace(ymin2,ymax2,len);

%% plot velocity 
figure("Name","Velocity")
plot3(plotVec*-0.8,span1,all_logs.plotOverLine_xm0_8.Velocity_0)
hold on; grid on;
plot3(plotVec*-0.4,span1,all_logs.plotOverLine_xm0_4.Velocity_0)
plot3(plotVec*0,span1,all_logs.plotOverLine_x0.Velocity_0)
plot3(plotVec*0.05,span2,all_logs.plotOverLine_x0_05.Velocity_0)
plot3(plotVec*0.1,span2,all_logs.plotOverLine_x0_1.Velocity_0)
plot3(plotVec*0.2,span2,all_logs.plotOverLine_x0_2.Velocity_0)
plot3(plotVec*0.3,span2,all_logs.plotOverLine_x0_3.Velocity_0)
plot3(plotVec*1,span2,all_logs.plotOverLine_x1.Velocity_0)
plot3(plotVec*1.4,span2,all_logs.plotOverLine_x1_4.Velocity_0)
legend('x=-0.8','x=-0.4','x=0','x=0.05','x=0.1','x=0.2','x=0.3','x=1','x=1.4')
title("Velocity 0 at different positions")
