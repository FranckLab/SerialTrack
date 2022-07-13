function plotF(xGrid,yGrid,zGrid,F_Grid_Vector)


subplot(3,3,1); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(1:9:end)); cb=colorbar; title('F11');
subplot(3,3,2); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(4:9:end)); cb=colorbar; title('F12');
subplot(3,3,3); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(7:9:end)); cb=colorbar; title('F13');
subplot(3,3,4); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(2:9:end)); cb=colorbar; title('F21');
subplot(3,3,5); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(5:9:end)); cb=colorbar; title('F22');
subplot(3,3,6); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(8:9:end)); cb=colorbar; title('F23');
subplot(3,3,7); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(3:9:end)); cb=colorbar; title('F31');
subplot(3,3,8); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(6:9:end)); cb=colorbar; title('F32');
subplot(3,3,9); 
scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),F_Grid_Vector(9:9:end)); cb=colorbar; title('F33');

end