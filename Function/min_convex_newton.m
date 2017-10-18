function [x,O_min] = min_convex_newton(O_fun, dd_f, d_f, A_eq,b_eq, O_neq,dd_neq,d_neq, x0)
%%
%input: O_fun,dd_f, d_f are functions 

%n_var = nargin(O_fun);
% parameter;
m = 100; M = 2e2; mu = 3; e = 10^(-5);
alpha = 0.4; beta = 0.6;
    
[len,n] = size(x0);
% if n ~= n_var
%     error('the number is not right!');
% end
b_len = length(b_eq);

x = x0;

while M/m > e
    while 1
        dd = dd_f(x)+dd_neq(x)/m;
        d = d_f(x) + d_neq(x)/m;
        A_dd = zeros(len*n);
        flag = 1;
        for i = 1:n
            for j = i:n
                if j == i
                    A_dd(1+(i-1)*len:i*len,1+(j-1)*len:j*len) = diag(dd(:,flag))/2;
                else
                    A_dd(1+(i-1)*len:i*len,1+(j-1)*len:j*len) = diag(dd(:,flag));
                end
                flag = flag + 1;
            end
        end
        A_dd = A_dd + A_dd';
        
        b = [-d; zeros(b_len,1)];
        A = [A_dd A_eq';
            A_eq zeros(b_len)];
        d_x_v = A\b;
        d_x = d_x_v(1:n*len);
        now_e = -(d'*d_x);
        if now_e < 0
            error('e<0?');
        end
        if now_e < e 
            break;
        end
        
        x_temp = reshape(x,n*len,1);
        temp = -x_temp./d_x;
        t_temp = 0.8*min(temp(logical(temp>0)));
        if length(t_temp) > 0
            t = min(t_temp,1);
        else
            t  =1;
        end
        
        if sum((x_temp + t*d_x)<0) > 0
            error('t is too big\n');
        end
        
        d_x_temp = reshape(d_x,len,n);
        O_F = O_fun(x)+O_neq(x)/m;
        O_F_test = O_fun(x+t*d_x_temp)+O_neq(x+t*d_x_temp)/m;
        count = 1;
        while O_F_test > (O_F-t*alpha*now_e)
            t = beta * t;
            O_F_test = O_fun(x+t*d_x_temp)+O_neq(x+t*d_x_temp)/m;
            count = count + 1;
            if count > 100
                break
            end
        end
        
        x = x+t*d_x_temp;
      
    end
    m = mu*m;
end

O_min = O_fun(x);
        