function adjust_fig(fig, ax, leg, x_label, y_label)

sizeTime = 1;

fig.Units = 'centimeters';
fig.Position(3:4) = sizeTime*[8, 5.25]; % 8cm x 5.25cm
fig.Color = [1.0, 1.0, 1.0]; % background color

% axes position and size: x y width height
ax.Units = 'centimeters';
ax.Position =  [1.4 + 0.5*sizeTime, 1, sizeTime*5.25,  4 ];
ax.LineWidth = 1;
ax.FontName = 'Times New Roman';
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';

ax.XLabel.String = x_label;
ax.XLabel.Units = 'normalized';
ax.XLabel.Position(1:2) =  [0.5, -0.125];
ax.XLabel.Interpreter = 'latex';
ax.XLim(1) = 0;

ax.YLabel.String = y_label;
ax.YLabel.Units = 'normalized';
ax.YLabel.Position(1:2) =  [-0.15, 0.5];
ax.YLabel.Interpreter = 'latex';
ax.YLim(1) = 0;

% tick
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength(1) = 0.02;

leg.Interpreter = 'latex';
leg.Location = 'best';
leg.Box = 'off';
end

%%
% fig = figure;
% ax = axes(fig);
% % blablabla
% 
% fig_new = figure;
% ax_new = copyobj(ax, fig_new);