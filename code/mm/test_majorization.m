% test for majorizing concave functions using quadratic surrogates. 

x = linspace(-1,1,2048);

epsilon = 10e-1;

%type = 'atan';
type = 'log2';
%type = 'log';
%type = 'abs';
%type = 'tukey';

% 1) function

g = 0;
switch type
    case 'log2'
        f = @(x)log(epsilon^2+x.^2) - log(epsilon^2);
        g = @(x)2*x./(epsilon^2 + x.^2);
    case 'log'
        f = @(x)log(epsilon+abs(x));% - log(epsilon);
    case 'atan'
        f = @(x)atan(abs(x)/epsilon);
    case 'abs'
        f = @(x)abs(x);
    case 'tukey'
        f = @(x) epsilon^2/6*(1-(1-(x/epsilon).^2).^3).*(abs(x)<=epsilon) ...
            + epsilon^2/6.*(abs(x)>epsilon);
    otherwise
        error(['Unknown : ' type]);
end
if g==0
    % derivative (using finite difference bc I am lazy
    g = @(x)(f(x+1e-9) - f(x))/1e-9;
end

% 2) quadratic majorizing function at point y
% h(x,y)  = f(y) + (x-y)*g(y) + alpha/2*(x-y)^2
%         = cst + (g(y) - alpha*y)*x + alpha/2 * x^2
% => rule to remove the linear part : alpha = g(y)/y

alpha = @(y)g(y)/y;
h = @(x,y)f(y) + (x-y)*g(y) + alpha(y)/2*(x-y).^2;

%
% 3) display for different values of y

Y = linspace(0.01, 0.2, 10);
clf; hold on;
plot(x,f(x), 'b');
for i=1:length(Y)
    y = Y(i); 
    t = (i-1)/(length(Y)-1);
    c = [t 1-t 0];
    plot(x,h(x,y)-f(x), 'color', c);
    plot(y,f(y), 'MarkerSize', 10, 'color', c);
end
axis([min(x) max(x) 0 max(f(x))*1.5]);