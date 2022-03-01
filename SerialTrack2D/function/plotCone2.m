function plotCone2(x,y,u,v)

p_sim_pw = [x(:),y(:),zeros(length(x(:)),1)];
u3_sim = [u(:),v(:),zeros(length(u(:)),1)];

% === Cone plot the solution ===
ps=p_sim_pw;  
u_mag = sqrt(u3_sim(:,1).^2 + u3_sim(:,2).^2 + u3_sim(:,3).^2); Umag_max = max(u_mag);
 %plot3(ps(:,1)+u(:,1),ps(:,2)+u(:,2),ps(:,3)+u(:,3),'k.'); hold on;
 
hc = coneplot(ps(:,1),ps(:,2),ps(:,3),u3_sim(:,1),u3_sim(:,2),u3_sim(:,3),0.04,'nointerp');  % default size: 0.04

try 
caxis([0,Umag_max]); fvc = repmat(u_mag(:)',[42 1]);
set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:));
hc.EdgeColor = 'none'; hc.AmbientStrength = 0.6; hc.DiffuseStrength = 0.75;  
hc.SpecularStrength = 0.4; hc.SpecularExponent = 3; 
h_color = colorbar; set(h_color, 'ylim', [0, Umag_max]); set(h_color,'fontsize',14);

lighting phong;

catch
end


title('Displacement (Unit: um)')

set(gca,'fontsize',14);  axis on;
ylabel('y (um)'); xlabel('x (um)');
grid minor; grid on; set(gcf,'color','w');

view(2);