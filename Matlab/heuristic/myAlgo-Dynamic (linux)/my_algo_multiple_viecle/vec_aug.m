function en_vec=vec_aug(vec,vec_d,n_add)
x_p=[vec(1:n_add+1);vec_d(1)];
x_pp=[vec(n_add+2:2*n_add+1);vec_d(2);vec(end)];
en_vec=[x_p;x_pp];