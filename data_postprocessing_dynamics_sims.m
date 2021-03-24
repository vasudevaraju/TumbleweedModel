function data_postprocessing_dynamics_sims(out,R,rc2,rr2)
%This is for chapter 2 of the thesis

y = out.simout.data;
t = out.simout.Time;
rm = out.simout1.Data;

[X,Y] = meshgrid(min(y(:,7))-3:1:max(y(:,7))+3,min(y(:,8))-3:1:max(y(:,8))+3);
Z = 0*X;%(4*exp(-((X-5).^2 + (Y).^2)./(6^2))-6*exp(-((X-5).^2 + (Y).^2)./(4^2)));%0*X;%cos(X./4) + 1;%-5*exp(-((X - 5).^2 + (Y - 0).^2)./(5^2));%(4*exp(-((X-5).^2 + (Y).^2)./(6^2))-6*exp(-((X-5).^2 + (Y).^2)./(4^2)));
figure;
surf(X,Y,Z,'FaceAlpha',0); hold on;
% [X,Y,Z] = sphere;
% surf(X*(rr2) + rc2(1), Y*(rr2) + rc2(2), Z*(rr2) + rc2(3) )
daspect([1,1,1]);

% cylR = 0.225;
% theta = [0:10*pi/180:360*pi/180].';
% z = [2.85:0.25:3.1].';
% [Theta,Z] = meshgrid(theta,z);
% CylR = cylR*ones(size(Theta,1),size(Theta,2));
% [X,Y,Z] = pol2cart(Theta,CylR,Z);
% surf(X,Y,Z,'FaceAlpha',0); hold on
% rm = zeros(size(t,1),18);

% [X, Y, Z] = sphere;
% surf(X*(0.2)+5, Y*(0.2), Z*0.2 - 0.1);
% surf(X*(0.5)+7, Y*(0.5), Z*0.5 - 0.1);
% rm = zeros(size(rm,1),18);
% th = [0:10*pi/180:360*pi/180]';
% r = 0.225;
% z_s = [2.8:0.1:3.1]';
% for i = 1:size(z_s,1)
% for j = 1:size(th,1)
% [X_s(i,j),Y_s(i,j),Z_s(i,j)] = pol2cart(th(j,1),r,z_s(i,1));
% end
% end
% surf(X_s,Y_s,Z_s,'FaceAlpha','0')
% hold on
daspect([1,1,1])
plot3(y(:,7),y(:,8),y(:,9),'LineWidth',1);%title('position of center of mass');
xlabel('x(m)');ylabel('y(m)');zlabel('z(m)'); hold on
draw_endpoint(y(end,:),t(end,:),rm(end,:), R); %xlim([min(y(:,7))-1,max(y(:,7))+1]); ylim([min(y(:,8))-1,max(y(:,8))+1]);xlabel('t(s)');xlim([min(t),max(t)])
figure;subplot(4,1,1);plot(t,y(:,10),'b','LineWidth',2);xlim([min(t),max(t)]);  ylabel(['$q_{_0}','$'],'interpreter','latex','FontSize',14); ylim([-1,1]);xlabel('t(s)');xlim([min(t),max(t)]);grid on; grid minor;
subplot(4,1,2);plot(t,y(:,11),'r','LineWidth',2);xlim([min(t),max(t)]);  ylabel(['$q_{_1}','$'],'interpreter','latex','FontSize',14); ylim([-1,1]);xlabel('t(s)');xlim([min(t),max(t)]);grid on; grid minor;
subplot(4,1,3);plot(t,y(:,12),'g','LineWidth',2);xlim([min(t),max(t)]);  ylabel(['$q_{_2}','$'],'interpreter','latex','FontSize',14);ylim([-1,1]);xlabel('t(s)');xlim([min(t),max(t)]);grid on; grid minor;
subplot(4,1,4);plot(t,y(:,13),'k','LineWidth',2);xlim([min(t),max(t)]);  ylabel(['$q_{_3}','$'],'interpreter','latex','FontSize',14); xlabel('t(s)'); ylim([-1,1]);xlim([min(t),max(t)]);grid on; grid minor;


figure;plot(t,y(:,1),t,y(:,2),t,y(:,3),'LineWidth',2);legend(['$\omega_{_x}','$'],['$\omega_{_y}','$'],['$\omega_{_z}','$'],'interpreter','latex','FontSize',14);xlabel('t(s)');ylabel(['$ \{^{^{\overline{O}}}\vec{\omega}^{^{\overline{T}}}\}_{_{\overline{T}}}\:(rad/s)','$'],'interpreter','latex','FontSize',14);xlim([min(t),max(t)]);grid on; grid minor;
figure;plot(t,y(:,4),t,y(:,5),t,y(:,6),'LineWidth',2);legend(['$\dot{x}_{_T}','$'],['$\dot{y}_{_T}','$'],['$\dot{z}_{_T}','$'],'interpreter','latex','FontSize',14);xlabel('t(s)');ylabel(['$\{{^{^{\overline{O}}}\vec{v}}_{_{T/O}}\}_{_{\overline{O}}}\:(m/s)','$'],'interpreter','latex','FontSize',14);xlim([min(t),max(t)]);grid on; grid minor;
figure;plot(t,y(:,7),t,y(:,8),t,y(:,9),'LineWidth',2);legend(['$x_{_T}','$'],['$y_{_T}','$'],['$z_{_T}','$'],'interpreter','latex','FontSize',14);xlabel('t(s)');ylabel(['$\{{\vec{r}}_{_{T/O}}\}_{_{\overline{O}}}\:(m)','$'],'interpreter','latex','FontSize',14);xlim([min(t),max(t)]);grid on; grid minor;

end

function draw_endpoint(y,t,rm, R)
control_cone = 0;
rm_dummy(:,1) =  rm(:,1);rm_dummy(:,2) =  rm(:,3);rm_dummy(:,3) =  rm(:,5);
rm_dummy(:,4) =  rm(:,2);rm_dummy(:,5) =  rm(:,4);rm_dummy(:,6) =  rm(:,6);
rm = rm_dummy;
r = R;
theta = linspace(0,2*pi,100);
phi = linspace(0,2*pi,100);
[Theta,Phi,R] = meshgrid(theta,phi,r) ;
[X,Y,Z] = sph2cart(Theta,Phi,R);
z_c = linspace(-r,r,50);
theta_c = linspace(0,2*pi,50);
[Theta_C, Z_C] = meshgrid(theta_c,z_c);
R_C = control_cone*Z_C;
[X_C,Z_C,Y_C] = pol2cart(Theta_C, R_C, Z_C);
i = 1;
    for j =1:size(X)
        for k = 1:size(X)
            pos_vec = rot_mat_OcT(y(i,:))*[X(j,k),Y(j,k),Z(j,k)]';
            X(j,k) = pos_vec(1);Y(j,k) = pos_vec(2);Z(j,k) = pos_vec(3);
        end
    end
    s1 = surf(X+y(i,7),Y+y(i,8),Z+y(i,9),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);hold on
    s2 = surf(X_C+y(i,7),Y_C+y(i,8),Z_C+y(i,9),'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
    OcT = rot_mat_OcT(y(i,:));
    l1 = line([y(i,7),y(i,7)+OcT(1,1)*(r)],[y(i,8),y(i,8)+OcT(2,1)*(r)],[y(i,9),y(i,9)+OcT(3,1)*(r)],'LineWidth',2,'Color','k');
    l2 = line([y(i,7),y(i,7)+OcT(1,1)*(-r)],[y(i,8),y(i,8)+OcT(2,1)*(-r)],[y(i,9),y(i,9)+OcT(3,1)*(-r)],'LineWidth',2,'Color','k');
    l3 = line([y(i,7),y(i,7)+OcT(1,2)*(r)],[y(i,8),y(i,8)+OcT(2,2)*(r)],[y(i,9),y(i,9)+OcT(3,2)*(r)],'LineWidth',2,'Color','k');
    l4 = line([y(i,7),y(i,7)+OcT(1,2)*(-r)],[y(i,8),y(i,8)+OcT(2,2)*(-r)],[y(i,9),y(i,9)+OcT(3,2)*(-r)],'LineWidth',2,'Color','k');
    l5 = line([y(i,7),y(i,7)+OcT(1,3)*(r)],[y(i,8),y(i,8)+OcT(2,3)*(r)],[y(i,9),y(i,9)+OcT(3,3)*(r)],'LineWidth',2,'Color','k');
    l6 = line([y(i,7),y(i,7)+OcT(1,3)*(-r)],[y(i,8),y(i,8)+OcT(2,3)*(-r)],[y(i,9),y(i,9)+OcT(3,3)*(-r)],'LineWidth',2,'Color','k');
    c1 = scatter3(y(i,7)+rm(i,1)*OcT(1,1),y(i,8)+rm(i,1)*OcT(2,1),y(i,9)+rm(i,1)*OcT(3,1),'LineWidth',1,'MarkerEdgeColor','g','MarkerFaceColor','g'); 
    c2 = scatter3(y(i,7)+rm(i,4)*OcT(1,1),y(i,8)+rm(i,4)*OcT(2,1),y(i,9)+rm(i,4)*OcT(3,1),'LineWidth',1,'MarkerEdgeColor','g','MarkerFaceColor','g');   
    c3 = scatter3(y(i,7)+rm(i,2)*OcT(1,2),y(i,8)+rm(i,2)*OcT(2,2),y(i,9)+rm(i,2)*OcT(3,2),'LineWidth',1,'MarkerEdgeColor','y','MarkerFaceColor','y');   
    c4 = scatter3(y(i,7)+rm(i,5)*OcT(1,2),y(i,8)+rm(i,5)*OcT(2,2),y(i,9)+rm(i,5)*OcT(3,2),'LineWidth',1,'MarkerEdgeColor','y','MarkerFaceColor','y');   
    c5 = scatter3(y(i,7)+rm(i,3)*OcT(1,3),y(i,8)+rm(i,3)*OcT(2,3),y(i,9)+rm(i,3)*OcT(3,3),'LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b');   
    c6 = scatter3(y(i,7)+rm(i,6)*OcT(1,3),y(i,8)+rm(i,6)*OcT(2,3),y(i,9)+rm(i,6)*OcT(3,3),'LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b');   
    center = scatter3(y(i,7),y(i,8),y(i,9),5,'MarkerEdgeColor','k','MarkerFaceColor','k');
    daspect([1,1,1]);
end

function OcT = rot_mat_OcT(y)
    TcO_val = [y(10)^2 + y(11)^2 - y(12)^2 - y(13)^2, 2*(y(11)*y(12) + y(10)*y(13)), 2*(y(11)*y(13) - y(10)*y(12));...
       2*(y(11)*y(12) - y(10)*y(13)), y(10)^2 - y(11)^2 + y(12)^2 - y(13)^2, 2*(y(12)*y(13) + y(10)*y(11));...
       2*(y(11)*y(13) + y(10)*y(12)), 2*(y(12)*y(13) - y(10)*y(11)), y(10)^2 - y(11)^2 - y(12)^2 + y(13)^2];
    OcT = transpose(TcO_val);
end