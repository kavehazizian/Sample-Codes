function en_vec=vec_aug_dyn(vec,vec_d,n_add)
x_p=[vec(1:n_add);vec_d(1)];
x_pp=[vec(n_add+1:2*n_add);vec_d(2);vec(end)];
en_vec=[x_p;x_pp];