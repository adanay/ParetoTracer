% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [] = mopvarlabels(psdim, sizes)
l = min(length(psdim), 3);
[n, nobj, na, naeq, nc, nceq] = val.valmopsizes(sizes, false);

for i = 1 : l
  if psdim(i) > n
    if psdim(i) > n + nobj
      if psdim(i) > n + nobj + na
        if psdim(i) > n + nobj + na + naeq
          if psdim(i) > n + nobj + na + naeq + nc
            if psdim(i) > n + nobj + na + naeq + nc + nceq
              if psdim(i) > n + nobj + na + naeq + nc + nceq + n % upper bounds
                label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$\varrho_{'', num2str(psdim(i) - n - nobj - na - naeq - nc - nceq - n), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex'''); 
              else % lower bounds
                label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$\rho_{'', num2str(psdim(i) - n - nobj - na - naeq - nc - nceq), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex''');
              end
            else % nonlin eqs
              label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$\lambda_{'', num2str(psdim(i) - n - nobj - na - nc), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex''');
            end
          else % nonlin ineqs
            label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$\gamma_{'', num2str(psdim(i) - n - nobj - naeq), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex''');
          end
        else % lin eqs
          label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$\lambda_{'', num2str(psdim(i) - n - nobj - na), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex''');
        end
      else % lin ineqs
        label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$\gamma_{'', num2str(psdim(i) - n - nobj), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex''');
      end
    else % objs
      label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$\alpha_{'', num2str(psdim(i) - n), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex''');
    end
  else % vars
    label(psdim, n, nobj, na, naeq, nc, nceq, i, '[''$x_{'', num2str(psdim(i)), ''}$''], ''fontsize'', 18, ''interpreter'', ''latex''');
  end
end
end

function [] = label(psdim, n, nobj, na, naeq, nc, nceq, i, str)
switch i
  case 1
    eval(['xlabel(' str ')']);
  case 2
    eval(['ylabel(' str ')']);
  case 3
    eval(['zlabel(' str ')']);
end
end

