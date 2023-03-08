function mod_dist=my_form(distal,n)
n_t=2*n+2;

fact=n_t./distal(:,2);
fact=fact.^3;
sd=(fact)-1;
mod_dist=(fact.*sd).*distal(:,3)+distal(:,3);
[df jh]=min(mod_dist);
distal(jh,:)