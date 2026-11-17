default(realprecision,38);
read("submission/search/mestre-own/mestre.gp");
b=[148,116,104,57,25,0];tt=82;qr=mesqr(concat(b+vector(6,i,tt),b-vector(6,i,tt)));print("OUR=",qr[2]);nr=(9*tt^2+211950)*x^4+(-2700*tt^2-63901710)*x^3+(-18*tt^4+396150*tt^2+6706476489)*x^2+(2700*tt^4-29575350*tt^2-284435346600)*x+9*tt^6-159200*tt^4+891699592*tt^2+4156297690000;print("RATIO=",qr[2]/nr);print("point residual=",((-224*tt^3-844*tt^2+900484*tt+2161725)/75)^2-subst(nr,x,(tt+703)/15));quit;
