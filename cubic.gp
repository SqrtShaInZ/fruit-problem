my_sqrt(n)=return(sqrtint(numerator(n))/sqrtint(denominator(n)));
gettors(E,n)=local(fulln,struct,P);[fulln,struct,P]=elltors(E);return(ellmul(E,P[1],fulln/n));
mkc(n)=return(ellinit([0,4*n^2+12*n-3,0,32*(n+3),0]));
mkc2(n)=return(ellinit([0,4*n^2+12*n-3,0,-128*(n+3),-128*(n+3)*(4*n^2+12*n-3)]));
mkc3(n)=return(ellinit([0,4*n^2+12*n-3,0,-320*n^2-1248*n-1104,-1024*n^4-7168*n^3-18944*n^2-24576*n-15040]));
mkc6(n)=return(ellinit([0,4*n^2+12*n-3,0,-640*n^3-6080*n^2-18528*n-18384,-2048*n^5-39936*n^4-281088*n^3-935936*n^2-1503744*n-941248]));
mkc_dispatch(n,ind)=if(ind==1,return(mkc(n)),if(ind==2,return(mkc2(n)),if(ind==3,return(mkc3(n)),return(mkc6(n)))));
exp1(n)=localprec(128);return(1/ellbsd(mkc(n))/2^ellrank(mkc(n))[3]);
exp2(n)=localprec(128);return(1/ellbsd(mkc2(n))/2^ellrank(mkc2(n))[3]);
exp3(n)=localprec(128);return(1/ellbsd(mkc3(n))/2^ellrank(mkc3(n))[3]);
exp6(n)=localprec(128);return(1/ellbsd(mkc6(n))/2^ellrank(mkc6(n))[3]);
exp_dispatch(n,ind)=localprec(128);if(ind==1,return(exp1(n)),if(ind==2,return(exp2(n)),if(ind==3,return(exp3(n)),return(exp6(n)))));
isogeny_point_aux(P)=return([P[1]/P[3]^2,P[2]/P[3]^3]);
\\ when all four curves have the same ellbsd, use this function
solve_isogeny_hard(n)=local(gen,x,y);gen=[findgen_worker(n,-600000,0)[1],[],[],[]];[x,y]=gen[1];\
gen[2]=ellsaturation(mkc2(n),[isogeny_point_aux(eval(ellisogeny(mkc(n),gettors(mkc(n),2))[2]))],1000)[1];\
gen[3]=ellsaturation(mkc3(n),[isogeny_point_aux(eval(ellisogeny(mkc(n),gettors(mkc(n),3))[2]))],1000)[1];\
gen[4]=ellsaturation(mkc6(n),[isogeny_point_aux(eval(ellisogeny(mkc(n),gettors(mkc(n),6))[2]))],1000)[1];\
comp2=ellheight(mkc(n),gen[1])<ellheight(mkc2(n),gen[2]);comp3=ellheight(mkc(n),gen[1])<ellheight(mkc3(n),gen[3]);\
if(comp2,if(comp3,return([n,6,gen[4]]),return([n,2,gen[2]])),if(comp3,return([n,3,gen[3]]),return([n,1,gen[1]])));
\\ when two curves have the same ellbsd, use this function
solve_isogeny(n,ind1,ind2)=if(ind1>ind2,my(t=ind1);ind1=ind2;ind2=t);local(P1,P2,E1,E2,isog);E1=mkc_dispatch(n,ind1);E2=mkc_dispatch(n,ind2);P1=ellrank(E1,3)[4];P2=ellrank(E2,3)[4];\
if(#P1==1&&#P2==1,if(ellheight(E1,P1[1])>ellheight(E2,P2[1]),return([n,ind1,P1[1]]),return([n,ind2,P2[1]])));\
if(#P1==0,if(#P2==0,return(solve_isogeny_hard(n)));if(ind2/ind1==3,return([n,ind1,map_III(n,E1,P2[1])]),return([n,ind1,map_II(n,E1,P2[1])])));\
[x,y]=P1[1];return([n,ind2,isogeny_point_aux(eval(isog=ellisogeny(E1,gettors(E1,ind2/ind1))[2]))]);
bestcurve_dispatch(n,ind)=local(x,y,e6,e3,e2,P2,P3,P6,ind_2isog,ind_3isog,e_orig,e_2isog,e_3isog);ind_2isog=if(ind%2==0,ind/2,ind*2);ind_3isog=if(ind%3==0,ind/3,ind*3);\
localprec(128);e_orig=exp_dispatch(n,ind);e_2isog=exp_dispatch(n,ind_2isog);e_3isog=exp_dispatch(n,ind_3isog);\
if(e_2isog-e_orig>2^-64&&e_3isog-e_orig>2^-64,return([n,ind,ellrank(mkc_dispatch(n,ind),3)[4][1]]));\
if(e_2isog-e_orig<=2^-64,if(e_3isog-e_orig<=2^-64,return(solve_isogeny_hard(n)),return(solve_isogeny(n,ind,ind_2isog))),return(solve_isogeny(n,ind,ind_3isog)));
bestcurve(n)=localprec(128);comp2=exp1(n)<exp2(n);comp3=exp1(n)<exp3(n);if(comp2,if(comp3,return(1),return(3)),if(comp3,return(2),return(6)));
gp_to_magma(res)=print(Str("ComputeGeneratorFull(", res[1], ", ", res[2], ", 2 : known_gens:=[[", res[3][1], ", ", res[3][2], "]]);"));
best(n)=if(ellrootno(mkc(n))==1,gp_to_magma(bestcurve_dispatch(n,bestcurve(n))),return(bestcurve(n)));
map_II(n,E,P)=[F,map]=ellisogeny(E,gettors(E,2));P=ellmul(ellinit(F),P,2);kill(x);Q=numerator(map[1]/map[3]-P[1]*map[3]);\
f=select(t->poldegree(t)==1,factor(Q)[,1])[1];P[1]=-polcoeff(f,0)/polcoeff(f,1);P[2]=my_sqrt(P[1]^3+E[2]*P[1]^2+E[4]*P[1]+E[5]);return(P);
map_III(n,E,P)=[F,map]=ellisogeny(E,gettors(E,3));P=ellmul(ellinit(F),P,3);kill(x);Q=numerator(map[1]-P[1]*map[3]^2);\
f=select(t->poldegree(t)==1,factor(Q)[,1])[1];P[1]=-polcoeff(f,0)/polcoeff(f,1);P[2]=my_sqrt(P[1]^3+E[2]*P[1]^2+E[4]*P[1]+E[5]);return(P);
map2(n,P)=return(map_II(n,mkc(n),P));
map3(n,P)=return(map_III(n,mkc(n),P));
map6(n,P)=return(map_II(n,mkc(n),map_III(n,mkc2(n),P)));
descend(E,limit,rankub)=local(poi,covers,Plist,hyperell,trans,i,P,x,y,hmat);poi=ellrank(E,abs(limit\200000))[4];return(poi);
\\ if(rankub==#poi,return(poi));covers=ell2cover(E);Plist=Vec([]);\
\\ [hyperell,trans]=covers[#covers];poi=hyperellratpoints(hyperell,abs(limit));for(i=1,#poi,[x,y]=poi[i];P=eval(trans);if(ellorder(E,P)==0,Plist=concat(Plist,[P]);\
\\ localprec(38);hmat=ellheightmatrix(E,Plist);if(abs(matdet(hmat))<1e-19,Plist=Plist[1..(#Plist)-1]));if(rankub==#Plist,return(Plist)));return(Plist);
heegner(n)=localprec(19);eh1=exp1(n);eh2=exp2(n);eh3=exp3(n);eh6=exp6(n);iso2=eh1>eh2;iso3=eh1>eh3;\
if(iso2,if(iso3,return(map6(n,ellheegner(mkc6(n)))),return(map2(n,ellheegner(mkc2(n))))),if(iso3,return(map3(n,ellheegner(mkc3(n)))),return(ellheegner(mkc(n)))));
try_isogeny_II(n,limit,rankub)=local(F,gens);F=mkc2(n);gens=descend(F,limit,rankub);if(gens==[],return([]));for(i=1,#gens,gens[i]=map_II(n,mkc(n),gens[i]));return(gens);
try_isogeny_III(n,limit,rankub)=local(F,gens);F=mkc3(n);gens=descend(F,limit,rankub);if(gens==[],return([]));for(i=1,#gens,gens[i]=map_III(n,mkc(n),gens[i]));return(gens);
try_isogeny_VI(n,limit,rankub)=local(F,gens);F=mkc6(n);gens=descend(F,limit,rankub);if(gens==[],return([]));for(i=1,#gens,gens[i]=map_III(n,mkc(n),map_II(n,mkc3(n),gens[i])));return(gens);
findgen_worker(n,limit=5*10^5,bound=1)=local(E,ranklb,rankub,Plist,gens,i,P,hmat);E=mkc(n);[ranklb,rankub]=ellrank(E)[1..2];if(bound&&ranklb!=rankub,local(ranklb2,rankub2);\
[ranklb2,rankub2]=ellrank(mkc2(n))[1..2];ranklb=max(ranklb,ranklb2);rankub=min(rankub,rankub2);print("rank bounded to [",ranklb,",",rankub,"]");\
if(ranklb!=rankub,localprec(19);if(ranklb<=1,if(rankub<4,rankub=ranklb+2*(abs(ellL1(E,ranklb))<10^-9),rankub=ellanalyticrank(E)[1]));print("Analytic rank is ",rankub)));\
localprec(38);if(rankub==0,return([]));Plist=[];gens=descend(E,limit,rankub);for(i=1,#gens,P=gens[i];Plist=concat(Plist,[P]);\
hmat=ellheightmatrix(E,Plist);if(abs(matdet(hmat))<1e-19,Plist=Plist[1..(#Plist)-1],print("Found point of height ",ellheight(E,P)," via direct search")));\
if(length(Plist)==rankub,return(Plist));gens=try_isogeny_II(n,limit,rankub);for(i=1,#gens,P=gens[i];Plist=concat(Plist,[P]);hmat=ellheightmatrix(E,Plist);\
if(#Plist>1&&abs(matdet(hmat))<1e-19,Plist=Plist[1..(#Plist)-1],print("Found point of height ",ellheight(E,P)," via 2-isogeny")));\
if(length(Plist)==rankub,return(Plist));gens=try_isogeny_III(n,limit,rankub);for(i=1,#gens,P=gens[i];Plist=concat(Plist,[P]);hmat=ellheightmatrix(E,Plist);\
if(#Plist>1&&abs(matdet(hmat))<1e-19,Plist=Plist[1..(#Plist)-1],print("Found point of height ",ellheight(E,P)," via 3-isogeny")));\
if(length(Plist)==rankub,return(Plist));gens=try_isogeny_VI(n,limit,rankub);for(i=1,#gens,P=gens[i];Plist=concat(Plist,[P]);hmat=ellheightmatrix(E,Plist);\
if(#Plist>1&&abs(matdet(hmat))<1e-19,Plist=Plist[1..(#Plist)-1],print("Found point of height ",ellheight(E,P)," via 6-isogeny")));\
if(length(Plist)==rankub,return(Plist),for(i=#Plist+1,rankub,Plist=concat(Plist,[[]]));if(rankub>=2,warning("Only partial basis has been found");return(Plist),\
if(limit<0,return([[]]),print1("[heegner] ");default(debug,1);return([heegner(n)]))));
findgen(n,reg_only=0,limit=5*10^5)=if(reg_only>=1,gens=ellrank(mkc(n),5);\
if(gens[2]==length(gens[4]),return(prod(i=1,gens[2],ellheight(mkc(n),gens[4][i]))),return(expectedheight(mkc(n),gens[2]))));return(findgen_worker(n,limit));
getxyz(N,P)=a=8*(N+3)-P[1]+P[2];b=8*(N+3)-P[1]-P[2];c=-8*(N+3)-2*(N+2)*P[1];g=gcd([a,b,c]);return([a/g,b/g,c/g]);
check_pos(n,P)=x=P[1];return(x<-4*(n+3)/(n+2)&&x^2+4*n*(n+3)*x+16*(n+3)^2>0);
solve_from_point(n,Preal,height_only=0)=if(Preal[1]>0,return([]));local(E,now_prec,P,tors,mul,P2,nowP);E=mkc(n);ellheight(E,Preal);now_prec=max(round(n*log(n)),default(realprecision));\
tors=gettors(E,6);mul=1;warning("Initial precision: ",now_prec);while(1,localprec(now_prec);P=Preal*1.0;nowP=P;P2=ellmul(E,P,2);mul=1;while(mul<now_prec*2/3,\
if(check_pos(n,nowP),break);if(check_pos(n,elladd(E,nowP,tors)),mul=-mul;break);mul=mul+2;nowP=elladd(E,nowP,P2);if(mul%100==99,print1(".")));if(abs(mul)>=now_prec*2/3,\
now_prec=round(now_prec*log(n));warning("recalculating with precision "now_prec),break));if(height_only<=-1,return(mul));add_tors=mul<0;mul=abs(mul);\
print1(mul,"*P");if(add_tors,print1("+tors"));localprec(38);h=3/2*ellheight(E,Preal)*mul^2-4*log(n)-9;print(" gives a positive solution, height>",h,", decimal digits>",h/log(10));\
if(height_only>=1,return(h));Preal=ellmul(E,Preal,mul);if(add_tors,return(getxyz(n,elladd(E,Preal,tors))),return(getxyz(n,Preal)));
find_j_range(i,h1,h2,pr,max_height)=my(a,b,c,d);a=h2;b=2*i*pr;c=i^2*h1-max_height;d=b^2-4*a*c;if(d<0,return([1,-1]));d=sqrt(d);return([ceil((-b-d)/(2*a)),floor((-b+d)/(2*a))]);
solve_from_two_egg_point(n,P1,P2,height_only=0)=my(E,mul1,mul2,h1,h2,pr,max_height,prec,i,j,cur_height,win_i,win_j,tors,P1R,iP1R,P2R);\
E=mkc(n);mul1=solve_from_point(n,P1,-1);mul2=solve_from_point(n,P2,-1);h1=ellheight(E,P1);h2=ellheight(E,P2);pr=ellheight(E,P1,P2);\
max_height=min(h1*mul1^2,h2*mul2^2);prec=min(abs(mul1),abs(mul2));localprec(prec);warning("precision: ",prec);P1R=P1*1.0;P2R=P2*1.0;win_i=0;win_j=0;tors=0;\
if(abs(max_height-h1*mul1^2)<10^-8,win_i=mul1;if(mul1<0,tors=1;win_i=-win_i;mul1=-mul1),win_j=mul2;if(mul2<0,tors=1;win_j=-win_j;mul2=-mul2));\
for(i=1,mul1-1,if(i==0,continue);iP1R=ellmul(E,P1R,i);jrange=find_j_range(i,h1,h2,pr,max_height);for(j=jrange[1],jrange[2],if(j==0||(i+j)%2==0,continue);\
cur_height=i^2*h1+j^2*h2+2*i*j*pr;target=elladd(E,iP1R,ellmul(E,P2R,j));if(check_pos(n,target),warning("New best point found with (i,j,tors)=(",i,",",j,",0)");\
max_height=cur_height;win_i=i;win_j=j;tors=0);target=elladd(E,target,[0,0]);if(check_pos(n,target),warning("New best point found with (i,j,tors)=(",i,",",j,",1)");\
max_height=cur_height;win_i=i;win_j=j;tors=1)));if(height_only>=1,return(win_i^2*h1+win_j^2*h2+2*win_i*win_j*pr));if(height_only<=-1,return([win_i,win_j,tors]));\
P=elladd(E,ellmul(E,P1,win_i),ellmul(E,P2,win_j));if(tors==1,return(getxyz(n,elladd(E,P,[0,0]))),return(getxyz(n,P)));
solve_from_egg_and_unbounded_point(n,P1,P2,height_only=0)=my(E,mul1,i,j,win_i,win_j,tors,P1R,P2R,max_height,cur_height);E=mkc(n);mul1=solve_from_point(n,P1,-1);\
win_i=abs(mul1);win_j=0;tors=mul1<0;h1=ellheight(E,P1);h2=ellheight(E,P2);pr=ellheight(E,P1,P2);max_height=h1*mul1^2;localprec(win_i);warning("precision: ",win_i);\
P1R=P1*1.0;P2R=P2*1.0;forstep(i=1,mul1,2,iP1R=ellmul(E,P1R,i);jrange=find_j_range(i,h1,h2,pr,max_height);for(j=jrange[1],jrange[2],if(j==0,continue);\
cur_height=i^2*h1+j^2*h2+2*i*j*pr;target=elladd(E,iP1R,ellmul(E,P2R,j));if(check_pos(n,target),warning("New best point found with (i,j,tors)=(",i,",",j,",0)");\
max_height=cur_height;win_i=i;win_j=j;tors=0);target=elladd(E,target,[0,0]);if(check_pos(n,target),warning("New best point found with (i,j,tors)=(",i,",",j,",1)");\
max_height=cur_height;win_i=i;win_j=j;tors=1)));if(height_only>=1,return(win_i^2*h1+win_j^2*h2+2*win_i*win_j*pr));if(height_only<=-1,return([win_i,win_j,tors]));\
P=elladd(E,ellmul(E,P1,win_i),ellmul(E,P2,win_j));if(tors==1,return(getxyz(n,elladd(E,P,[0,0]))),return(getxyz(n,P)));
solve_from_two_point(n,P1,P2,height_only=0)=if(P1[1]>0,if(P2[1]<0,return(solve_from_egg_and_unbounded_point(n,P2,P1,height_only)),return([])),\
if(P2[1]<0,return(solve_from_two_egg_point(n,P1,P2,height_only)),return(solve_from_egg_and_unbounded_point(n,P1,P2,height_only))));
solvexyz(n,height_only=0)=local(E,gens);if(n%2==1,return([]));E=mkc(n);gens=findgen(n);if(#gens==0,return([]));if(#gens>2,error("Not yet implemented"));\
if(#gens==2,gens=ellsaturation(E,gens,1000);return(solve_from_two_point(n,gens[1],gens[2],height_only)));return(solve_from_point(n,gens[1],height_only));
prepare_db(lb=1,ub=10000,limit=-5*10^5)=local(i,gens);for(i=lb,ub,print("finding generator for ",i);gens=findgen_worker(i,limit);write("cubic_db.txt",i," ",gens));
prepare_db_list(lst,limit=-5*10^5)=local(i,gens,n);for(i=1,#lst,n=lst[i];print("finding generator for ",n);gens=findgen_worker(n,limit);print("Generators of E_n: ",gens));
\\ grep -E -x '[0-9]+ \[\[\]\]' cubic_db.txt
\\
\\ SetSeed(1);
\\ SetClassGroupBounds("GRH");
\\ E := EllipticCurve([0, 147161149, 0, -142934720408016, -16846971914709437693120]);
\\ HyperE := TwoDescent(E : RemoveTorsion := true)[1];
\\ ///// below are 6-descent codes
\\ Crv3, mapA := ThreeDescentByIsogeny(E);Crv, map6to3 := SixDescent(HyperE, Crv3[1]);P6 := PointSearch(Crv, 10^6 : OnlyOne := true)[1];P3_internal := map6to3(P6);comps := Components(mapA[1]);P3_fixed := Domain(comps[1]) ! Eltseq(P3_internal);P_intermediate := comps[1](P3_fixed);P_intermediate_fixed := Domain(comps[2]) ! Eltseq(P_intermediate);PE := comps[2](P_intermediate_fixed);P := Saturation([PE], 1000 : TorsionFree := true)[1];P in E;CanonicalHeight(P);
\\ ///// below are 4-descent codes
\\ //Crv4 := FourDescent(HyperE)[1];Ps := PointsQI(Crv4, 10^7 : OnlyOne := true)[1];A, mapA := AssociatedEllipticCurve(Crv4 : E := E);P := Saturation([mapA(Ps)], 1000 : TorsionFree := true)[1];P in E;CanonicalHeight(P);
\\ ///// below are 12-descent codes
\\ //Crv4 := FourDescent(HyperE)[1];Crv3 := ThreeDescentByIsogeny(E)[1];Crvs12, maps12 := TwelveDescent(Crv3, Crv4);Crv12 := Crvs12[2]; map12 := maps12[2];P12 := PointSearch(Crv12, 10^6 : OnlyOne := true)[1];P4 := map12(P12);A, mapA := AssociatedEllipticCurve(Crv4 : E := E);P := Saturation([mapA(P4)], 1000 : TorsionFree := true)[1];P in E;CanonicalHeight(P);
