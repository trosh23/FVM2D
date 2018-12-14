x = linspace(0,1,10);
y = linspace(0,2,10);

[xq,yq] = meshgrid(0:0.01:1,0:0.02:2);

vq1 = griddata(x,y,data,xq,yq);
vq2 = griddata(x,y,data2,xq,yq);

hold on
mesh(xq,yq,vq1);
mesh(xq,yq,vq2);
hold off