%this file creates miscellaneous commands that are used in the NR multivariate method

%create string of dx inputs
params='x(1)';
if length(x)>1
  for i=2:length(x)
    params=strcat(params,',x(',num2str(i),')');
  end
end

%create string of fprintf initial input
msglog='Iteration %d: x1=%.6f';
if length(x)>1
  for i=2:length(x)
    msglog=strcat(msglog,', x',num2str(i),'=%.6f');
  end
end

msglog=strcat(msglog,', dx1=%.6f');
if length(x)>1
  for i=2:length(x)
    msglog=strcat(msglog,', dx',num2str(i),'=%.6f');
  end
end

msglog=strcat(msglog,', err1=%.6f');
if length(x)>1
  for i=2:length(x)
    msglog=strcat(msglog,', err',num2str(i),'=%.6f');
  end
end

msglog=strcat(msglog,'\n');