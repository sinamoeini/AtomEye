% compare_angular.m %

load angular.out;

clf;
plot(angular(:,1), angular(:,2), angular(91,1), angular(91,2), 'ro');

grid on;
xlabel('\theta [degree]');
ylabel('A(\theta) [eV/A]');

print -depsc compare_angular.eps;
