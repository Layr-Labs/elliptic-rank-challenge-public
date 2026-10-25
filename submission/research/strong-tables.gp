{
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([7*t^2+(309),0,(-91)*(t^2-1)*(t^2-25)*(t^2-121),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),9999));Einf=ellinit([Mod(7,p),0,Mod(-91,p),0,0]);apinf=if(#Einf,p+1-ellcard(Einf),9999);print(1," ",p," ",concat(vap,[apinf])));
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([13*t^2+(3711),0,(-1189)*(t^2-25)*(t^2-121)*(t^2-361),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),9999));Einf=ellinit([Mod(13,p),0,Mod(-1189,p),0,0]);apinf=if(#Einf,p+1-ellcard(Einf),9999);print(2," ",p," ",concat(vap,[apinf])));
}
quit;
