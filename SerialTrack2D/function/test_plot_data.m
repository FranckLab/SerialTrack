
close all; clear; clc

%% Demo I
x=1:1:10; 

y1=x+10; y2=x+20; y3=x+30; y4=x+40; y5=x+50; y6=x+60; y7=x+70; y8=x+80; y9=x+90;

datax{1}=x; datax{2}=x; datax{3}=x; datax{4}=x; 
% datax{5}=x; datax{6}=x; datax{7}=x; datax{8}=x; datax{9}=x;

datay{1}=y1; datay{2}=y2; datay{3}=y3; datay{4}=y4; 
% datay{5}=y5; datay{6}=y6; datay{7}=y7; datay{8}=y8; datay{9}=y9;

display_names=['1';'2';'3';'4'];

[fig, ax, lgd] = plot_data(datax,datay,display_names);

adjust_fig(fig, ax, lgd, '','');

axis([0,11,0,55]);


%% Demo II
tempx = 1:0.1:10;
tempy = sin(tempx);

fig=figure; ax=axes;
plot(tempx,tempy,'o-');
lgd = []; %lgd = legend('I am legend','position','northeastOutside');
xlabel='$x$ axis';
ylabel='$y$ axis';
adjust_fig(fig, ax, lgd, xlabel,ylabel);
 
axis([0,11,-1,1])

