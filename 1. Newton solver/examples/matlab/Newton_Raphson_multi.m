function x=Newton_Raphson_multi(x0,f1,epsilon,maxiter=100,print=true,J_inv=[])
  %inputs: x0 is initial x-value, f1 is the symbolic function, epsilon is the precision
  if isempty(J_inv)
    J=jacobian(f1);
    J_inv=inv(J);
  end  
  dx=-J_inv*f1;
  dx=matlabFunction(dx);
  
  x=x0;
  iter=0;
  
  NR_multi_commands;
  
  %loop
  while eval(strcat('sum(abs(dx(',params,')))'))>epsilon*length(x0)
    if print==true
      NR_multi_printing;
    end  
    x=x+eval(strcat('dx(',params,')')); %update x=x+dx(x)
    iter=iter+1;
    if iter>maxiter
      h=msgbox('Max Iterations Reached');
      break
    end  
  end
  %final print
  if print==true
    NR_multi_printing;
  end  
end  