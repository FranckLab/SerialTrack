function plotScatter3(x,y,z,u,v,w)

x=x(:); y=y(:); z=z(:); u=u(:); v=v(:); w=w(:);

subplot(2,2,1); 
scatter3(x,y,z,ones(length(x),1),u); cb=colorbar; title('Displacement x');
subplot(2,2,2); 
scatter3(x,y,z,ones(length(x),1),v); cb=colorbar; title('Displacement y');
subplot(2,2,3); 
scatter3(x,y,z,ones(length(x),1),w); cb=colorbar; title('Displacement z');
subplot(2,2,4); 
u_mag = sqrt(u.^2 + v.^2 + w.^2);
scatter3(x,y,z,ones(length(x),1),u_mag); cb=colorbar; title('|Disp|');


