function c = heat(nu,n,dt,nsteps,nplot,mflag)

//
// c = heat(nu,n,dt,nsteps,mflag)
//
// Compute the solution of the heat equation u_t - nu u_xx = f(x)
// with homogeneous Dirichlet boundary conditions, starting from the
// initial solution  u0(x). 
//
// nu  is the diffusion coefficient.
// n   is the number of internal nodes on the interval [0,1].
// dt  is the time step.
// nsteps  is the number of time steps.
// nplot   plot the solution every nplot time steps.
//
// mflag = 0  forward Euler time-stepping
// mflag = 1  backward Euler time-stepping
// mflag = 2  Crank-Nicholson time-stepping
//

  h = 1/(n+1);    // space step

// computation of the transfer matrix
  
  bet = nu*dt/h^2;
  D = sparse([1:n;1:n]',2*bet*ones(n,1),[n,n]);
  L = sparse([2:n;1:n-1]',-bet*ones(n-1,1),[n,n]);
  U = sparse([1:n-1;2:n]',-bet*ones(n-1,1),[n,n]);
  A = L+D+U;
  
// initial solution
  x = h:h:1-h;
  u = u0(x)';
  xx = 0:h:1;
  yy = u0(xx);
  
  select mflag

  case 0 then
  // Forward Euler

  case 1 then
  // Backward Euler
  M = speye(n,n)+A;

  case 2 then
  // Cranck-Nicholson
  M = speye(n,n)+0.5*A;

  end;

  for iter = 1:nsteps
  
    if (modulo(iter,nplot)==0)
      yy(2:n+1) = u(1:n,1)';
      plot(xx,yy,'r',x,f(x),rect=[0,-0.5,1,1.5])
      set(gca(),"auto_clear","on")
    end;
    
    select mflag
    
    case 0 then
    // Forward Euler
    u = u - A*u + dt*f(x)';

    case 1 then
    // Backward Euler
    u = M\(u + dt*f(x)');
    
    case 2 then
    // Cranck-Nicholson
    u = M\(u - 0.5*A*u + dt*f(x)');
    
    end;
    
  end;

  c = u;

endfunction

//The subfunctions defininf the source term and intial conditions

 function f=f(x)
 //f=0.0;
 f=exp(-200*(x-0.25).^2);
 endfunction
 
 function u0=u0(x)
 //u0=4.*x.*(1-x);
 //u0=20.*max(0,(x-0.25).*(0.75-x));
 //u0=sin(19*%pi*x);
 //u0=sin(%pi*x)+0.5*sin(4*%pi*x);
 u0=0.0*x;
 endfunction
