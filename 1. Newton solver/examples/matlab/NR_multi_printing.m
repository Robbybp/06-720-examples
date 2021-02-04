%this file prints the iterations for the NR multivariate method
Dx=eval(strcat('dx(',params,')'));
aDx=abs(Dx);

params1=num2str(x(1));
if length(x)>1
  for i=2:length(x)
    xi=num2str(x(i));
    params1=strcat(params1,', ',xi);
  end
end

params2=num2str(Dx(1));
if length(x)>1
  for i=2:length(x)
    dxi=num2str(Dx(i));
    params2=strcat(params2,', ',dxi);
  end
end

params3=num2str(aDx(1));
if length(x)>1
  for i=2:length(x)
    abs_dx=num2str(aDx(i));
    params3=strcat(params3,', ',abs_dx);
  end
end

msgprint=strcat('fprintf(''',msglog,''', ',num2str(iter),', ',params1,', ',params2,', ',params3,')');
eval(msgprint);
fflush(stdout); 