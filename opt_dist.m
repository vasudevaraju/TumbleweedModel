function [point_on_surf, smallest_dist] = opt_dist(r_TO)
u0 = [r_TO(1);r_TO(2)];
[point_on_surf, smallest_dist] = fminunc(@(u)surf_dist(u,r_TO),u0);
end

function z_c = surf_zcoord(x_c,y_c, r_TO)
% if r_TO(3) >= 0
%     z_c = sqrt(x_c^2 - 2*y_c);
% else
%     z_c = -sqrt(x_c^2 - 2*y_c);
% end
z_c = 0;%(4*exp(-((x_c-5)^2 + (y_c)^2)/(6^2))-6*exp(-((x_c-5)^2 + (y_c)^2)/(4^2)));%(-1*exp(-((x_c - 10)^2 + (y_c - 10)^2)/2^2) -5*exp(-((x_c - 10)^2 + (y_c - 15)^2)/5^2) + 5*exp(-((x_c - 10)^2 + (y_c + 10)^2)/5^2) - 3*exp(-((x_c - 20)^2 + (y_c +5)^2)/4^2) - 5*exp(-((x_c - 40)^2 + (y_c -10)^2)/6^2));
end

function dist_val = surf_dist(u,r_TO)
x = u(1); y = u(2);
% if r_TO(2) >= 0
%     y = sqrt(0.225^2 - x^2);
% else
%    y = -sqrt(0.225^2 - x^2);
% end
z = surf_zcoord(x,y, r_TO);
% z = r_TO(3);
dist_val = norm([x;y;z] - r_TO);
end