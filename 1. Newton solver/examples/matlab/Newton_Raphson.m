function x=Newton_Raphson(x0,f1,epsilon,maxiter=100)
  %f1=0
  df1=diff(f1);
  df2=diff(df1);
  f1=matlabFunction(f1);
  df1=matlabFunction(df1);
  df2=matlabFunction(df2);
  convergence=abs(f1(x0)*df2(x0)/df1(x0)^2)
  x=x0;
  iter=0;
  Error=1;
  while Error>epsilon && iter <= maxiter 
     dxi=-f1(x)/df1(x);
     Error=abs(dxi);
     fprintf('Iteration %d: xi=%.6f, dxi=%.6f, err=%.6e\n', iter, x, dxi, Error);
     x_old=x;
     x=x_old+dxi;
     iter=iter+1;
  end
end
