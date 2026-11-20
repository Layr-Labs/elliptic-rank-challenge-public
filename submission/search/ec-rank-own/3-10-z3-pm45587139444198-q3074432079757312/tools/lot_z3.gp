\\ lot_z3.gp (build-lottery, 19 Nov 2026): Z/3, y^2 + a x y + b y = x^3 (point (0,0) of order 3); (a,b) ~ (u a, u^3 b), so the curve is
\\ E_t: y^2 + x y + t y = x^3 with t = b / a^3.  Sieve lines read 'z3 p q' with p = b and q = a^3 (so t = p/q); integral model
\\ [q, 0, p q^2, 0, 0]; ellminimalmodel(ellinit(famcurve(p/q))) is the minimal model of y^2 + a x y + b y = x^3.
famcurve(t) = [1, 0, t, 0, 0];
