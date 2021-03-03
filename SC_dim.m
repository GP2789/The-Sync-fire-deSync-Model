function [ ] = SC_dim( N1, N2, N3 )
%SC_DIM Summary of this function goes here
%   Detailed explanation goes here
U1 = 1:N1+N2+N3;
U2 = [N1 N1*(1:N2+N3)];
U3 = [N1*N2 N1*N2*(1:N3)];

T1 = 1:N1+N2+N3;
T2 = N1:N1+N2+N3;
T3 = N1+N2:N1+N2+N3;

figure(); hold on; box on;
plot(T1,U1,'k-','linewidth',3)
plot(T2,U2,'k--','linewidth',3)
plot(T3,U3,'k:','linewidth',3)

xticks([0 N1 N1+N2 N1+N2+N3])
xlim([1 N1+N2+N3])
yticks([N1 N1*N2 N1*N2*N3])
ylim([0 max(U3)])
set(gca, 'fontsize', 20, 'fontname','arial');
xlabel('unique nodes')
ylabel('unique time points')

end

