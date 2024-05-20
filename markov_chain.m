% 5 years matrix
z=[1	1.97	3.89	3.24	6.69	12.61	137.36	1349.46];
p1=[0.874871	0.119119	0.004004	0	0	0	0.002006	0
0.06	0.877996	0.06	0	0	0	0.002004	0
0.004004	0.119119	0.874871	0	0	0	0.002006	0
0	0	0	0.874871	0.119119	0.004004	0.002006	0
0	0	0	0.06	0.877996	0.06	0.002004	0
0	0	0	0.004004	0.119119	0.874871	0.002006	0
0.021	0.021	0.021	0.021	0.021	0.021	0.85	0.024
0	0	0	0	0	0	0.242	0.758];

row_sums = sum(p1, 2);  % Calculate the sum of each row
% initial distibution
pini=[0.044	0.412	0.044	0.044	0.412	0.044	0	0];
% working age distibution

sum_mat = zeros(size(p1)); % initialize sum matrix
for i = 1:9
    sum_mat = sum_mat + p1^i;
end
%population share	0.097	0.303	0.097	0.097	0.303	0.097	0.0076	0.0003
pop=pini*1/9*sum_mat

income = pop * z';

sh_i=pop .* z/income*100;

pop_a=[pop(1:2), sum(pop(3:4)),pop(5:end)]'*100
sh_ia=[sh_i(1:2), sum(sh_i(3:4)),sh_i(5:end)]'

[pop'*100, sh_i']
% Data

data=[-10 2;-40 14;60 14;40 34; 10 20; 1 13; 0.1 6];
sum(data)
disp('Data vs Model') % 20+13+6
[data, pop_a, sh_ia]

% Annualizing the matrix
pa=zeros(size(p1));
pa(1:6,1:6)=p1(1:6,1:6)^(1/5);
pa(7:8,7:8)=p1(7:8,7:8)^(1/5);
pa(1:6,7)=mean(p1(1:6,7))/5;%1-mean(sum(pa(1:6,1:6), 2));
pa(7,1:6)=mean(p1(7,1:6))/5;

row_sums = sum(pa, 2);  % Calculate the sum of each row

pa= pa ./ row_sums;
row_sums = sum(pa, 2)


sum_mat = zeros(size(pa)); % initialize sum matrix
for i = 1:45
    sum_mat = sum_mat + pa^i;
end
%population share	0.097	0.303	0.097	0.097	0.303	0.097	0.0076	0.0003
pop=pini*1/45*sum_mat;

income = pop * z';

sh_i=pop .* z/income*100;

pop_a=[pop(1:2), sum(pop(3:4)),pop(5:end)]'*100
sh_ia=[sh_i(1:2), sum(sh_i(3:4)),sh_i(5:end)]'

[pop'*100, sh_i']
% Data

disp('Data vs Model') % 20+13+6
[data, pop_a, sh_ia]











p = [0.8740 0.1190 0.0040; 0.0600 0.8780 0.0600; 0.0040 0.1190 0.8740];

row_sums = sum(p, 2);  % Calculate the sum of each row
pr = p ./ row_sums;  % Divide each element of p by the corresponding row sum

pr.^100;


pini=[0.125803	0.248395	0.125803	0.125803	0.248395	0.125803	0	0
];
pini*p1^5

pini=[0.1236    0.2441+0.08    0.1236    0.1236    0.2441+0.08    0.1236 	0	0
];
pini*p1^50



p1^10000