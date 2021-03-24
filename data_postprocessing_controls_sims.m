function data_postprocessing_controls_sims(out,R)
y = out.simout.data;
t = out.simout.Time;
rm = out.simout1.Data;
x1_comp = out.simout6.Data;
x3_comp = out.simout9.Data;
yref = out.simout2.data;
x1_track = out.simout14.Data;
x3_track = out.simout15.Data;
log = out.logsout{2}.Values.Data;

[X,Y] = meshgrid(min(y(:,7))-3:1:max(y(:,7))+3,min(y(:,8))-3:1:max(y(:,8))+3);
Z = 0*X;% (4*exp(-((X-5).^2 + (Y).^2)./(6^2))-6*exp(-((X-5).^2 + (Y).^2)./(4^2)));%0*X;%cos(X./4) + 1;%-5*exp(-((X - 5).^2 + (Y - 0).^2)./(5^2));%(4*exp(-((X-5).^2 + (Y).^2)./(6^2))-6*exp(-((X-5).^2 + (Y).^2)./(4^2)));
figure;
% surf(X,Y,Z,'FaceAlpha',0); hold on;

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
p2 = plot3(yref(:,1),yref(:,2),yref(:,3),'LineWidth',1);view(0,90);hold on
p1 = plot3(y(:,7),y(:,8),y(:,9),'LineWidth',1);
hold on;%title('position of center of mass'); 

grid on; grid minor
daspect([1,1,1]);

xlabel('x(m)');ylabel('y(m)');zlabel('z(m)'); hold on
draw_endpoint(y(end,:),t(end,:),rm(end,:), R); xlim([min(y(:,7))-3,max(y(:,7))+3]); ylim([min(y(:,8))-3,max(y(:,8))+3]);
legend([p1(1),p2(1)],['$\underline{x}_{_1}','$'],['$\underline{x}_{_{1,c}}','$'],'interpreter','latex','FontSize',14);

figure;subplot(4,1,1);hold on; plot(t,log(:,10));plot(t,y(:,10),'y','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$q_{_0}','$'],'interpreter','latex','FontSize',14); ylim([-1,1]);grid on; grid minor
subplot(4,1,2);hold on; plot(t,log(:,11));plot(t,y(:,11),'r','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$q_{_1}','$'],'interpreter','latex','FontSize',14); ylim([-1,1]);grid on; grid minor
subplot(4,1,3);hold on; plot(t,log(:,12));plot(t,y(:,12),'g','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$q_{_2}','$'],'interpreter','latex','FontSize',14);ylim([-1,1]);grid on; grid minor
subplot(4,1,4);hold on; plot(t,log(:,13));plot(t,y(:,13),'k','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$q_{_3}','$'],'interpreter','latex','FontSize',14); xlabel('t(s)'); ylim([-1,1]);grid on; grid minor

figure;subplot(3,1,1);hold on; plot(t,log(:,1));plot(t,y(:,1),'y','LineWidth',2);xlim([min(t),max(t)]);xlabel('t(s)'); ylabel(['$\omega_{_x}\:(rad/s)','$'],'interpreter','latex','FontSize',14); ylim([-1,1]);grid on; grid minor
subplot(3,1,2);hold on; plot(t,log(:,2));plot(t,y(:,2),'r','LineWidth',2);xlim([min(t),max(t)]);xlabel('t(s)'); ylabel(['$\omega_{_y}\:(rad/s)','$'],'interpreter','latex','FontSize',14); ylim([-1,1]);grid on; grid minor
subplot(3,1,3);hold on; plot(t,log(:,3));plot(t,y(:,3),'g','LineWidth',2);xlim([min(t),max(t)]);xlabel('t(s)'); ylabel(['$\omega_{_z}\:(rad/s)','$'],'interpreter','latex','FontSize',14);ylim([-1,1]);grid on; grid minor

figure;hold on; plot(t,log(:,1),t,log(:,2),t,log(:,3));plot(t,y(:,1),t,y(:,2),t,y(:,3),'LineWidth',2);xlabel('t(s)');legend(['$\omega_{_x}','$'],['$\omega_{_y}','$'],['$\omega_{_z}','$'],'interpreter','latex','FontSize',14);xlabel('t(s)');ylabel(['${\underline{x}}_{_3}\:(rad/s)','$'],'interpreter','latex','FontSize',14);xlim([min(t),max(t)]);
grid on; grid minor

figure;plot(t,log(:,7),t,log(:,8),t,log(:,9)-0.5);plot(t,y(:,7),t,y(:,8),t,y(:,9)-0.5,'LineWidth',2);legend(['$x','$'],['$y','$'],['$z','$'],'interpreter','latex','FontSize',14);xlabel('t(s)');ylabel(['${\underline{x}}_{_1}\:(m)','$'],'interpreter','latex','FontSize',14);xlim([min(t),max(t)]);
grid on; grid minor

figure;subplot(3,1,1);plot(t,rm(:,1),'b','LineWidth',2);xlim([min(t),max(t)]);ylabel(['$l_{_1}\:(m)','$'],'interpreter','latex','FontSize',14); xlabel(['$t\:(s)','$'],'interpreter','latex','FontSize',14);ylim([min(rm(:,1))-0.1,max(rm(:,1))+0.1]);grid on; grid minor
subplot(3,1,2);plot(t,rm(:,3),'r','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$l_{_2}\:(m)','$'],'interpreter','latex','FontSize',14);xlabel(['$t\:(s)','$'],'interpreter','latex','FontSize',14); ylim([min(rm(:,3))-0.1,max(rm(:,3))+0.1]);grid on; grid minor
subplot(3,1,3);plot(t,rm(:,5),'g','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$l_{_3}\:(m)','$'],'interpreter','latex','FontSize',14);xlabel(['$t\:(s)','$'],'interpreter','latex','FontSize',14);ylim([min(rm(:,5))-0.1,max(rm(:,5))+0.1]);grid on; grid minor

figure;subplot(3,1,1);plot(t,rm(:,7),'b','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$\dot{l}_{_1}\:(m/s)','$'],'interpreter','latex','FontSize',14);xlabel(['$t\:(s)','$'],'interpreter','latex','FontSize',14); ylim([min(rm(:,7))-0.1,max(rm(:,7))+0.1]);grid on; grid minor
subplot(3,1,2);plot(t,rm(:,9),'r','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$\dot{l}_{_2}\:(m/s)','$'],'interpreter','latex','FontSize',14);xlabel(['$t\:(s)','$'],'interpreter','latex','FontSize',14); ylim([min(rm(:,9))-0.1,max(rm(:,9))+0.1]);grid on; grid minor
subplot(3,1,3);plot(t,rm(:,11),'g','LineWidth',2);xlim([min(t),max(t)]); ylabel(['$\dot{l}_{_3}\:(m/s)','$'],'interpreter','latex','FontSize',14);xlabel(['$t\:(s)','$'],'interpreter','latex','FontSize',14);ylim([min(rm(:,11))-0.1,max(rm(:,11))+0.1]);grid on; grid minor

figure; plot(t,vecnorm(x1_comp,2,2),'LineWidth',2);xlabel('t(s)');ylabel(['$|\:\overline{\underline{x}}_{_1}\:|\:(m)','$'],'interpreter','latex','FontSize',14);grid on; grid minor; hold on

figure; plot(t,vecnorm(x3_comp,2,2),'LineWidth',2);xlabel('t(s)');ylabel(['$|\:\overline{\underline{x}}_{_3}\:|\:(rad/s)','$'],'interpreter','latex','FontSize',14);grid on; grid minor

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