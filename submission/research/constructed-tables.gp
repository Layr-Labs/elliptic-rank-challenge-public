{
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([7*t^2+(309),0,(-91)*(t^2-1)*(t^2-25)*(t^2-121),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),0));print(1," ",p," ",vap));
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([49*t^2+(-1521),0,(4160)*(t^2-1)*(t^2-49)*(t^2-169),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),0));print(2," ",p," ",vap));
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([8*t^2+(276),0,(-294)*(t^2-1)*(t^2-25)*(t^2-81),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),0));print(3," ",p," ",vap));
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([58*t^2+(-1444),0,(6156)*(t^2-1)*(t^2-16)*(t^2-361),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),0));print(4," ",p," ",vap));
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([3*t^2+(46813),0,(-2772224)*(t^2-1)*(t^2-169)*(t^2-225),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),0));print(5," ",p," ",vap));
forprime(p=2,1000,vap=vector(p);for(j=0,p-1,t=Mod(j,p);E=ellinit([53*t^2+(4588),0,(-143325)*(t^2-1)*(t^2-64)*(t^2-144),0,0]);vap[j+1]=if(#E,p+1-ellcard(E),0));print(6," ",p," ",vap));
}
quit;
