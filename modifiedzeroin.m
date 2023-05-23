function [root,info] = modifiedzeroin3032054936(func,Int,params)
% On input:
%   func is a function handle
%   Int is an interval, [Int.a, Int.b] containing a root
%   params is an object containing params.root_tol, params.func_tol and
%   params.maxit
% The algorithm will terminate once the interval containing the root is at
% most params.root_tol in length or the function value at the current
% iterate is at most params.func_tol in absolute value.
%
% The assumption of this method is that f(a)*f(b) < 0 and there is exactly
% one sign change on the interval [a,b], i.e. exactly one root in [a,b]
%
% On output:
%   root is the computed root
%   info has info.flag which is 0 on success and 1 otherwise
%
%
% initialize params
% set info.flg to indicate failure if root is not found
info.flg = 1;
% set N = max number of iterations
N = params.maxit;
% set r_tol to interval length tolerance
r_tol = params.root_tol;
% set f_tol to function value tolerance
f_tol = params.func_tol;

% initialize interval [a,b] and midpoint c
a = Int.a;
b = Int.b;
c = (a+b)/2;
% compute initial function values
fa = func(a);
fb = func(b);
fc = func(c);
% set initial x_1,x_2,x_3    
x(1) = a;
x(2) = b;
x(3) = c;
% set initial function values fx_1, fx_2, fx_3
fun(1) = fa;
fun(2) = fb;
fun(3) = fc;
% set gflag = 0 to indicate starting with x_1,x_2,x_3 from bisection
gflag = 0;
% set iflag = 0 to indicate number of consecutive IQI iterations
iflag = 0;
root = x(3);

% iterate until max number of iterations reached or root found
for k = 1:N
% check if we have found the root
    if (abs(b-a)<= r_tol) || (abs(fun(3)) <= f_tol)
% x_3 will always be our most current root approximation, return it and
% set info.flg = 0 indicating success
        root = x(3);
        info.flg = 0;
        return;
    end
% when gflag > 0, indicates IQI failed and 3 bisection steps are taken
    if gflag > 0
        for i = 1:3
% check if previously computed midpoint, c, is outside of new interval
% set during consecutive IQI steps 
            if (c<a) || (c>b)
                c = (a+b)/2;
                fc = func(c);
            end

 % find where sign changes, and set new a,b,c           
            if (sign(fc)*sign(fb)>0)

                b = c;
                fb = fc;
            else
                a = c;
                fa = fc;
            end
            c = (a+b)/2;
            fc = func(c);

            
% check if root is found during bisection steps
            if (abs(b-a)<= r_tol) || (abs(fc)<= f_tol)
                root = c;
                info.flag = 0;
                return;
            end
        end
 % reset x_1,x_2,x_3 after bisection steps
        x(1) = a;
        x(2) = b;
        x(3) = c;
% reset fx_1,fx_2,fx_3
        fun(1) = fa;
        fun(2) = fb;
        fun(3) = fc;
% reset gflag = 0 and iflag = 0, indicating x_1,x_2,x_3 are not IQI 
% iterates, and starting with x_1,x_2,x_3 from bisection steps
        gflag = 0;
        iflag = 0;
    end

% set x_4 = IQI(x_1,x_2,x_3)
    iq_1 = (((fun(2)*fun(3))/((fun(1)-fun(2))*(fun(1)-fun(3))))*x(1));
    iq_2 = (((fun(1)*fun(3))/((fun(2)-fun(1))*(fun(2)-fun(3))))*x(2));
    iq_3 = (((fun(1)*fun(2))/((fun(3)-fun(1))*(fun(3)-fun(2))))*x(3));

    x(4) =  iq_1 + iq_2 + iq_3;        

% check if x_4 is in interval containing root [a,b]
    if (x(4) < a || x(4) > b)
% when x_4 is not in [a,b], set gflag = 1, indicating IQI failed, then
% break and do 3 bisection steps at top of loop and start over
        gflag = 1;
        continue
    end
% we know that x_4 is in [a,b], compute its value

    fun(4) = func(x(4));

% set x_1=x_2, x_2=x_3, x_3=x_4
    x(1) = x(2); fun(1) = fun(2);
    x(2) = x(3); fun(2) = fun(3);
    x(3) = x(4); fun(3) = fun(4);
% set new endpoints of interval, first need to find relative order of
% points x_1, x_2, x_3. Then see where sign changes.
%
% variable rm represents the right most point of the new x_1,x_2,x_3 with
% frm being the function value of the right most point
%
% variable lm represents the left most point of the new x_1,x_2,x_3 with
% flm being the function value of the left most point
    if (x(3) > x(1)) && (x(3) > x(2))
        rm = x(3); frm = fun(3);
        if (x(2)>x(1))
            lm = x(1); flm = fun(1);
        else
            lm = x(2); flm = fun(2);
        end
    end
    if (x(2) > x(1)) && (x(2) > x(3))
        rm = x(2); frm = fun(2);
        if (x(3)>x(1))
            lm = x(1); flm = fun(1);
        else
            lm = x(3); flm = fun(3);
        end
    end
    if (x(1)>x(2)) && (x(1) > x(3))
        rm = x(1); frm = fun(1);
        if (x(3)>x(2))
            lm = x(2); flm = fun(2);
        else
            lm = x(3); flm = fun(3);
        end
    end
    
 % At this point we have our right most point as rm and left most point as
 % lm. We also know that x_1,x_2,x_3 are contained in the interval [lm,rm]
 % which is contained in [a,b].
 % Now we set new interval that contains x_1,x_2,x_3 and is contained
 % in [a,b], based on where sign changes.
    if ((sign(frm)*sign(fb))>0)
 % When the function value of the right most point has the same sign as the
 % function value of previous b, this indicates that the root must be
 % between the previous a and the right most point, so set b_new = rm.
        b=rm; fb=frm;

        if ((sign(flm)*sign(fa))>0)
 % When the function value of the left most point has the same sign as the
 % function value of the previous a, and we know from previous if statement
 % that the right most point had the same sign as the previous b, then we
 % know that the left most point and right most point must have different
 % signs, thus the root is between them, so set a_new = lm.
            a=lm; fa=flm;
        end
 % When the function value of the left most point has a different sign as 
 % the function value of the previous a, then the root must be between a
 % and the left most point, but the new interval must contain x_1,x_2,x_3
 % so we leave a unchanged, i.e. set a_new = a.
    else
        a=lm; fa=flm;
 % When the function value of the right most point has a different sign as
 % the function value of the previous b, the root must be between the right
 % most point and b, but the new interval must contain x_1,x_2,x_3 so we
 % set a_new = lm and leave b unchanged, i.e. set b_new = b.
    end

% If abs(f(x_3)) has not decreased by factor of 2 within 3 consecutive IQI 
% iterations, when iflag > 1 indicating that x_1,x_2,x_3 are all 
% IQI iterates break and do 3 bisection steps at top of loop, start over
    if ((2*abs(fun(3))) >= (abs(fun(1)))) && (iflag > 1)
        gflag = 1;
        continue
    end

% add 1 to iflag to keep track of number of consecutive IQI iterations
    iflag = iflag + 1;


end