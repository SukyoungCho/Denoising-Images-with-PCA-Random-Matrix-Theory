%% Normal distribution
STD = 1;
MEAN = 0;
x = -5:0.1:5;
f = ( 1/(STD*sqrt(2*pi)) ) * exp(-0.5*((x-MEAN)/STD).^2 );

plot (x,f);
title('Normal distribution');

%% SY - Normal
x = -5:0.1:5;
mu = 0; % 평균
sigma = 1; % 표준편차

y = normpdf(x, mu, sigma); % pdf for normal dist

plot(x, y);
xlabel('x'); 
ylabel('y');
title('Normal distribution');
%% Normal with hist
%Generate a sample of size 100 from a normal distribution 
%with mean 0 and variance 1.
r = normrnd(0,1,100,1);
hist(r);
%title('Normal distribution');
%histfit uses fitdist to fit a distribution to data. 
%Use fitdist to obtain parameters used in fitting.
%pd = fitdist(r,'Normal');

%% T distribution

%STD = 1;
%MEAN = 0;
nu= 5; % df
%f = ( 1/(STD*sqrt(2*pi)) ) * exp(-0.5*((x-MEAN)/STD).^2 );
%t = (1/sqrt(pi*v))*(factorial((v+1)/2-1)/factorial(v/2-1))*(1+t^2/v)^(-(v+1)/2);

t = tpdf(x,nu);

plot (x,t);
%xlabel('x'); 
%ylabel('t');
%title('T distribution');
histogram(t,15);

%y = tpdf(x,nu) returns the probability density function (pdf)
%of the Student's t distribution at each of the values in x using
%the corresponding degrees of freedom in nu. x and nu can be vectors,
%matrices, or multidimensional arrays that have the same size.
%A scalar input is expanded to a constant array with the same dimensions
%as the other inputs.

%% Hist - T distribution

r = trnd(10*ones(1,100));
hist(r);


%% T dist - upgrade
x = [-5:.1:5];
y1 = tpdf(x,5);   % For nu = 5
%y2 = tpdf(x,25);  % For nu = 25
%y3 = tpdf(x,50);  % For nu = 50

figure;
plot(x,y1,'Color','black','LineStyle','-')
%hold on
%plot(x,y2,'Color','red','LineStyle','-.')
%plot(x,y3,'Color','blue','LineStyle','--')
%legend({'nu = 5','nu = 25','nu = 50'})
%hold off
title('T distribution (nu = 5)');

%% Chi-squared
chi = chi2pdf(x,10); % I assigned v as 4. it could be any degree of freedom.
plot (x,chi);
xlabel('x'); 
ylabel('Chi');
title('Chi-Square distribution df = 10');


%% Chi-square hist
V = 10;                                             % Create Data
R = chi2rnd(V, 1, 100);                             % Create Data
X = histcounts(R,20);                               % Create Data
RNCF = @(v) norm(X - chi2pdf((1:length(X)),v));     % Residual Norm Cost Function
Ve = fminsearch(RNCF, rand);                        % Extimate Parameter
figure(1)
bar((1:length(X)), X/sum(X))

%% Chi-square curve
V = 10;                                             % Create Data
R = chi2rnd(V, 1, 100);                             % Create Data
X = histcounts(R,20);                               % Create Data
RNCF = @(v) norm(X - chi2pdf((1:length(X)),v));     % Residual Norm Cost Function
Ve = fminsearch(RNCF, rand);                        % Extimate Parameter
figure(1)
plot((1:length(X)), chi2pdf((1:length(X)),Ve), '-r')
title('Chi-Square Distribution');

%% Random Matrix and Eigen value

% 5*5
pd = makedist('Normal');
r = random(pd,[5,5]);
A = 1/25*r*r'; %1/n*E*E'
e = eig(A);
histogram(e,5);
title('5 X 5 random matrix');

% 10*10
pd = makedist('Normal');
r = random(pd,[10,10]);
A = 1/(10*10)*r*r'; %1/n*E*E'
e = eig(A);
histogram(e,5);
title('10 X 10 random matrix');

% 50*50
pd = makedist('Normal');
r = random(pd,[50,50]);
A = 1/(50*50)*r*r'; %1/n*E*E'
e = eig(A);
histogram(e,6);
title('50 X 50 random matrix');

% 250*250
pd = makedist('Normal');
r = random(pd,[250,250]);
A = 1/(250*250)*r*r'; %1/n*E*E'
e = eig(A);
histogram(e,10);
title('250 X 250 random matrix');

