function N=shapeFn2d(gp)
% zi and eta are column vectors
% zi=gp(:,1);
% eta=gp(:,2);
N=[(1/4)*(1-gp(:,1)).*(1-gp(:,2))     (1/4)*(1+gp(:,1)).*(1-gp(:,2)) ...
                    (1/4)*(1+gp(:,1)).*(1+gp(:,2))          (1/4)*(1-gp(:,1)).*(1+gp(:,2))];