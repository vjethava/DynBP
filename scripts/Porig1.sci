// plot the b vs p for node 0x0 at N1xN2_sY where Y is the seed
b00 = read('B00.txt', -1, 3);
p00 = read('P00.txt', -1, 1);
f = scf(1);
plot(b00(:, 1), p00 , 'r-d', b00(:, 1), b00(:, 2), 'g--.',b00(:, 1), b00(:, 3), 'b--*');
legend(['TrueP' 'ApproxP' 'RgnP'],4);
xtitle("Evolution of marginal probability at node (0,0) with time", "time" , "P(x00 = 0)");
xs2eps(1, 'b_vs_p.eps');
// plot D||P
DvP = read('DVP.txt', -1, 2);
g = scf(2);
plot(DvP(:, 1) , DvP(:, 2), 'b--d');
xtitle("Evolution of KL-divergence between true probability and DynBP approximation", "time", "KL(b||p)");
xs2eps(2, 'KL.eps');
h=scf(3); 
C = read('BbyP.txt', -1, 4); 
plot(C(:, 1), C(:, 2),'b-o' , C(:,1), C(:,3), 'r-x', C(:,1),  C(:, 4),'g-*');
xtitle("ratio of difference between beliefs and true probability normalized w."+...
       " r.t. true probability", "time", "(b-p)/p")
legend("min", "avg", "max", 4); 
xs2eps(3, 'ratio.eps')
