% plotangular.m %

load angular.out;

clf;
plot(angular(:,1), angular(:,2));

grid on;
xlabel('\theta [degree]');
ylabel('A(\theta) [eV/A]');

print -depsc angular.eps;
