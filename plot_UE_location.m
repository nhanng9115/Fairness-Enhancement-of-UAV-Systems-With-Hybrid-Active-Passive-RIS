function plot_UE_location(sys,V,RIS_vec,n_chan)

u = sys.u; r = sys.r; K = sys.K; dmin = sys.dmin; dmax = sys.dmax;

figure
plot(r(1,1),r(2,1),'gs','MarkerFaceColor','c','MarkerSize',40); hold on;

for k = 1:K
    plot(u(1,:),u(2,:),'rp','MarkerFaceColor','r','MarkerSize',12,'LineWidth',2); hold on;
end

if nargin >= 2
    linestyle = {'ko','g^','cs','md','bh','bs','--m'};

    for rr = 1:length(RIS_vec)
        RIS = RIS_vec(rr);
        plot(V(1,1:n_chan,rr),V(2,1:n_chan,rr),linestyle{RIS+1},'MarkerSize',6,'LineWidth',2); hold on;
    end
end
xlim([dmin dmax])
ylim([dmin dmax])
xlabel('x [m]')
ylabel('y [m]')