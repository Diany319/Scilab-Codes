function u=secondorderlinear(L,R,LV,RV,m,graphflag)

 // A GENERAL LINEAR SECOND ORDER ODE SOLVER
 // a(x)u"(x)+ b(x)u'(x)+c(x)u(x)=f(x)
 // on the interval L<=x<=R
 // with two dirichlet boundary conditions 
 // u(L)=LV,  u(R)=RV
 // for m internal nodes
 // for any functions a(x), b(x), c(x), f(x)
 // specified below as subfunctions
 // Solve directly the linear system AU=F
 //
 // L=x0, R=x(m+1), LV=u(x0), RV=u(x(m+1))
 // m=number of internal gridpoints
 // if graphflag>0 graph the solution, otherwise don't.
 
 h=(R-L)/(m+1);
 
 for i=1:m
     xi=R+i*h;
     // diagonal elements of A
     d(i)=h^2*c(xi)-2*a(xi);
 end
 
 for i=1:m-1
     xi=R+i*h;
     // subdiagonal elements
     ld(i)=a(xi)-h*b(xi)/2; 
     // superdiagonal elements
     ud(i)=a(xi)+h*b(xi)/2;
 end
 D=sparse([1:m;1:m]',d,[m,m]);
 UDIAG=sparse([1:(m-1);2:m]',ud,[m,m]);
 LDIAG=sparse([2:m;1:(m-1)]',ld,[m,m]);
 //The A matrix
 A=LDIAG+D+UDIAG;
 
 //The rhs F
 F=ones(m,1);
 for i=2:m-1
     xi=R+i*h;
     F(i)=h^2*f(xi);
 end
 xi=R+h;
 F(1)=h^2*f(xi)-(a(xi)-h*b(xi)/2)*LV;
 xi=L-h;
 F(m)=h^2*f(xi)-(a(xi)+h*b(xi)/2)*RV;
 
 u=A\F; //Solution
 
 if graphflag>0
     xx=L:(R-L)/(m+1):R;
     yy=ones(m+2);
     yy(1)=LV;
     yy(2:m+1)=u';
     yy(m+2)=RV;
     plot(xx,yy)
 end
 
 endfunction
 
 
 //The subfunctions define the functions a,b,c,f
 
 function a=a(x)
 a=-0.01;
 endfunction

 function b=b(x)
 b=1.0;
 endfunction

 function c=c(x)
 c=0.0;
 endfunction

 function f=f(x)
 f=0.0;
 endfunction

