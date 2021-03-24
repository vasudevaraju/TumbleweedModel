
function d3_animate(t,y,rm,control_cone,R,rock_cents,rock_radi,rock_there)
% input arguments : time t - n x 1 vector (should be evenly spaced)
% states y - n x 13 vector 
% where y(i,:) = [omegax,omegay,omegaz,u_TO,v_TO,w_TO,x,y,z,q0,q1,q2,q3]
% mass position vector rm - n x 6 vector
% where rm(i,:) = [x_m1,y_m3,z_m5,x_m2,y_m4,z_m6] 
% control_cone - scalar (just set to some positive number or 0)
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
%cylinder_render(y);
surf_render(y);
rm_dummy(:,1) =  rm(:,1);rm_dummy(:,2) =  rm(:,3);rm_dummy(:,3) =  rm(:,5);
rm_dummy(:,4) =  rm(:,2);rm_dummy(:,5) =  rm(:,4);rm_dummy(:,6) =  rm(:,6);
rm = rm_dummy;
r = R;
rock_rendered1 = false;
rock_rendered2 = false;
count = 1;
%surf_render(y);
theta = linspace(0,2*pi,100);
phi = linspace(0,2*pi,100);
[Theta,Phi,r1] = meshgrid(theta,phi,r) ;
[X,Y,Z] = sph2cart(Theta,Phi,r1);
z_c = linspace(-R,R,50);
theta_c = linspace(0,2*pi,50);
[Theta_C, Z_C] = meshgrid(theta_c,z_c);
R_C = control_cone*Z_C;
[X_C,Z_C,Y_C] = pol2cart(Theta_C, R_C, Z_C);
% v = VideoWriter('square_track_noise','MPEG-4');
% open(v)
for i = 1:2:size(t,1)
    for j =1:size(X)
        for k = 1:size(X)
            pos_vec = rot_mat_OcT(y(i,:))*[X(j,k),Y(j,k),Z(j,k)]';
            X(j,k) = pos_vec(1);Y(j,k) = pos_vec(2);Z(j,k) = pos_vec(3);
        end
    end
%     s1 = surf(X(1:25,1:25)+y(i,7),Y(1:25,1:25)+y(i,8),Z(1:25,1:25)+y(i,9),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);hold on
%     s2 = surf(X(1:25,26:50)+y(i,7),Y(1:25,26:50)+y(i,8),Z(1:25,26:50)+y(i,9),'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
%     s3 = surf(X(1:25,51:75)+y(i,7),Y(1:25,51:75)+y(i,8),Z(1:25,51:75)+y(i,9),'EdgeColor','none','FaceColor','b','FaceAlpha',0.5);
%     s4 = surf(X(1:25,76:100)+y(i,7),Y(1:25,76:100)+y(i,8),Z(1:25,76:100)+y(i,9),'EdgeColor','none','FaceColor','k','FaceAlpha',0.5);
%     s5 = surf(X(51:75,1:25)+y(i,7),Y(51:75,1:25)+y(i,8),Z(51:75,1:25)+y(i,9),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);hold on
%     s6 = surf(X(51:75,26:50)+y(i,7),Y(51:75,26:50)+y(i,8),Z(51:75,26:50)+y(i,9),'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
%     s7 = surf(X(51:75,51:75)+y(i,7),Y(51:75,51:75)+y(i,8),Z(51:75,51:75)+y(i,9),'EdgeColor','none','FaceColor','b','FaceAlpha',0.5);
%     s8 = surf(X(51:75,76:100)+y(i,7),Y(51:75,76:100)+y(i,8),Z(51:75,76:100)+y(i,9),'EdgeColor','none','FaceColor','k','FaceAlpha',0.5)
if rock_there
if y(i,7) < count*50 - 20
    if ~rock_rendered1
        clf('reset');hold on;
        xlim([y(i,7)-10,y(i,7)+10]);
        ylim([y(i,8)-10,y(i,8)+10]);
        zlim([0,2*R+2]);
        surf_render(y);
       rock_render(rock_cents(:,:,count),rock_radi(count,:));
        rock_rendered1 = true;
    end
elseif count*50 - 20 < y(i,7) && y(i,7) < count*50 + 20
      if ~rock_rendered2 && (count + 1)<=size(rock_radi,1)
        rock_render(rock_cents(:,:,count:count+1),rock_radi(count:count+1,:));
        rock_rendered2 = true;
      elseif ~rock_rendered2 && (count + 1)>size(rock_radi,1)
          rock_render(rock_cents(:,:,count),rock_radi(count,:));
          rock_rendered2 = true;
      end
else
    rock_rendered1 = false;
    rock_rendered2 = false;
    count = count + 1; 
end
end
%     view([0,-0.25,0.2]);
   %view([0,-1,0.5]);
%    view([0,30]);
view([45,45]);
    s1 = surf(X+y(i,7),Y+y(i,8),Z+y(i,9),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);hold on
    s2 = surf(X_C+y(i,7),Y_C+y(i,8),Z_C+y(i,9),'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
    OcT = rot_mat_OcT(y(i,:));
    l1 = line([y(i,7),y(i,7)+OcT(1,1)*(R)],[y(i,8),y(i,8)+OcT(2,1)*(R)],[y(i,9),y(i,9)+OcT(3,1)*(R)],'LineWidth',2,'Color','k');
    l2 = line([y(i,7),y(i,7)+OcT(1,1)*(-R)],[y(i,8),y(i,8)+OcT(2,1)*(-R)],[y(i,9),y(i,9)+OcT(3,1)*(-R)],'LineWidth',2,'Color','k');
    l3 = line([y(i,7),y(i,7)+OcT(1,2)*(R)],[y(i,8),y(i,8)+OcT(2,2)*(R)],[y(i,9),y(i,9)+OcT(3,2)*(R)],'LineWidth',2,'Color','k');
    l4 = line([y(i,7),y(i,7)+OcT(1,2)*(-R)],[y(i,8),y(i,8)+OcT(2,2)*(-R)],[y(i,9),y(i,9)+OcT(3,2)*(-R)],'LineWidth',2,'Color','k');
    l5 = line([y(i,7),y(i,7)+OcT(1,3)*(R)],[y(i,8),y(i,8)+OcT(2,3)*(R)],[y(i,9),y(i,9)+OcT(3,3)*(R)],'LineWidth',2,'Color','k');
    l6 = line([y(i,7),y(i,7)+OcT(1,3)*(-R)],[y(i,8),y(i,8)+OcT(2,3)*(-R)],[y(i,9),y(i,9)+OcT(3,3)*(-R)],'LineWidth',2,'Color','k');
    c1 = scatter3(y(i,7)+rm(i,1)*OcT(1,1),y(i,8)+rm(i,1)*OcT(2,1),y(i,9)+rm(i,1)*OcT(3,1),'LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g'); 
    c2 = scatter3(y(i,7)+rm(i,4)*OcT(1,1),y(i,8)+rm(i,4)*OcT(2,1),y(i,9)+rm(i,4)*OcT(3,1),'LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b');   
    c3 = scatter3(y(i,7)+rm(i,2)*OcT(1,2),y(i,8)+rm(i,2)*OcT(2,2),y(i,9)+rm(i,2)*OcT(3,2),'LineWidth',2,'MarkerEdgeColor','y','MarkerFaceColor','y');   
    c4 = scatter3(y(i,7)+rm(i,5)*OcT(1,2),y(i,8)+rm(i,5)*OcT(2,2),y(i,9)+rm(i,5)*OcT(3,2),'LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b');   
    c5 = scatter3(y(i,7)+rm(i,3)*OcT(1,3),y(i,8)+rm(i,3)*OcT(2,3),y(i,9)+rm(i,3)*OcT(3,3),'LineWidth',2,'MarkerEdgeColor','w','MarkerFaceColor','w');   
    c6 = scatter3(y(i,7)+rm(i,6)*OcT(1,3),y(i,8)+rm(i,6)*OcT(2,3),y(i,9)+rm(i,6)*OcT(3,3),'LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b');   
    center = scatter3(y(i,7),y(i,8),y(i,9),5,'MarkerEdgeColor','k','MarkerFaceColor','k');

    %traj_point = trajectory(y(i,7));
    %desired = scatter3(traj_point(1),traj_point(2),0.5,5,'MarkerEdgeColor','r','MarkerFaceColor','r');
    %axis([-1,10,-2,2,-4,4]);
    daspect([1,1,1]);
%       xlim([-2;12]);ylim([-2;12]);
    drawnow;
  
    xlim([y(i,7)-5,y(i,7)+5]);
    ylim([y(i,8)-5,y(i,8)+5]);
    %zlim([y(i,9)-5,y(i,9)+3]);
    zlim([0,2*R+2]);
        xlabel('metres');
    ylabel('metres');
    titletext = num2str(t(i));
    title("Simulation time : " + titletext + "s");
       xlabel('metres');
    ylabel('metres');
%     frame = getframe(gcf);
%     writeVideo(v,frame);
    if i ~= size(t,1)
        delete(s1);
        delete(s2);
%         delete(s3);
%         delete(s4);
%         delete(s5);
%         delete(s6);
%         delete(s7);
%         delete(s8);
        delete(l1);
        delete(l2);
        delete(l3);
        delete(l4);
        delete(l5);
        delete(l6);       
        delete(c1);
        delete(c2);
        delete(c3);
        delete(c4);
        delete(c5);
        delete(c6); 
    end
end
%  close(v)
end

function OcT = rot_mat_OcT(y)
    TcO_val = [y(10)^2 + y(11)^2 - y(12)^2 - y(13)^2, 2*(y(11)*y(12) + y(10)*y(13)), 2*(y(11)*y(13) - y(10)*y(12));...
       2*(y(11)*y(12) - y(10)*y(13)), y(10)^2 - y(11)^2 + y(12)^2 - y(13)^2, 2*(y(12)*y(13) + y(10)*y(11));...
       2*(y(11)*y(13) + y(10)*y(12)), 2*(y(12)*y(13) - y(10)*y(11)), y(10)^2 - y(11)^2 - y(12)^2 + y(13)^2];
    OcT = transpose(TcO_val);
end

function traj_point = trajectory(x)
T = 0  ; 
if x<T
        traj_point(1) = x;
        traj_point(2) = 0;
elseif x>=T
        traj_point(1) = x;
        traj_point(2) = sin(0.2*x);
end
end

function rock_render(rock_cents,rock_radi)
% persistent limit
% persistent count
% if isempty(limit)
%     limit = 50;
%     count = 1;
% end
% if (y(7) + 20 < limit)
 rcs = [];   
[X,Y,Z] = sphere;
for i = 1:size(rock_radi,1)
    for j = 1:size(rock_radi,2)
    rcs = [rcs;surf(X*(rock_radi(i,j)) + rock_cents(1,j,i), Y*(rock_radi(i,j)) + rock_cents(2,j,i), Z*(rock_radi(i,j)) + rock_cents(3,j,i))]; 
    end
end
end

function surf_render(y)
xs = min(y(:,7))-5:max(y(:,7))+5;
ys = min(y(:,8))-5:max(y(:,8))+5;
[Xs,Ys] = meshgrid(xs,ys);
Zs = zeros(size(Xs,1),size(Xs,2));
surf(Xs,Ys,Zs,'FaceAlpha',0);hold on;
% xs = min(y(:,7))-10:max(y(:,7))+10;
% ys = min(y(:,8))-10:max(y(:,8))+10;
% [Xs,Ys] = meshgrid(xs,ys);
% Zs = (4*exp(-((Xs-5).^2 + (Ys).^2)./(6^2))-6*exp(-((Xs-5).^2 + (Ys).^2)./(4^2)));
%  surf(Xs,Ys,Zs,'FaceAlpha',0);hold on;
end

function cylinder_render(y)
cylR = 0.225;
theta = [0:10*pi/180:360*pi/180].';
z = [2.85:0.25:3.1].';
[Theta,Z] = meshgrid(theta,z);
R = cylR*ones(size(Theta,1),size(Theta,2));
[X,Y,Z] = pol2cart(Theta,R,Z);
surf(X,Y,Z,'FaceAlpha',0)
end
    