function [point_on_surf, smallest_dist] = dist_calc(r_TO)
 coder.extrinsic('opt_dist');
[point_on_surf, smallest_dist] = opt_dist(r_TO);
end