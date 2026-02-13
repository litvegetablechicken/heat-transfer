clear
close all
[lzl,lzsa,lza,lzm,lzv,Torc,Tex,h,Tw_d,Tw_do] = datareduction('garcia_parallel_flow1123coorelation3-1');
%[lzl,lzsa,lza,lzm,lzv,Torc,Tex,h] = datareductionsimple('garcia_counter_flowcorrelation3');
target_value = 5.56;
fluid = 'R22';
din = 7.9E-3;
differences = abs(Torc(:, 1) - target_value);
[~, closest_index] = min(differences);
closest_value = Torc(closest_index, 1);
closest_corresponding_value = Torc(closest_index, 2);
A = pi * din * closest_value;
T1 = closest_corresponding_value+273.15;
P1 = 11.5*100;
Tcon = 25+273.15;
Tsc = 0;
G = 289;
cross_area = pi*((din/2.0)^2.0);
m =  G *cross_area;
Ets = 0.8;
[Wt, Wp, Qeva, Qcon, H2, P2, H4, P4, H1] = ORC(fluid, T1, P1, Tcon, Tsc, Ets, m);
cost = economic(Wt, Wp, A*2);
data = cost;
% Convert structure to double array
fields = fieldnames(data);  % Get field names
values = zeros(1, numel(fields));  % Pre-allocate array

for i = 1:numel(fields)
    values(i) = data.(fields{i});  % Assign values from structure to array
end

disp(values);  % Display resulting array
