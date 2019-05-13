A = [ 3/10 3/5 0 ; 1/2 0 1; 4/10 4/5 0];
W = [ 3/5 0 ; 0 1 ; 4/5 0];

H = W\A

At = transpose(A);
Ht = transpose(H);

Wt = Ht\At

W = transpose(Wt)