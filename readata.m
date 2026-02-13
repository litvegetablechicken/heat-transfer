function unique_data = readata(excel,sheet1,sheet3,sheet5,lzl,lza,lzv,varargin)
if length(varargin) == 4
    sheet2 = varargin{1};
    sheet4 = varargin{2};
    lzsa =  varargin{3};
    lzm = varargin{4};
    Tzsaex = read(excel,sheet2);
    Tzmex = read(excel,sheet4); 
    Tlzsuex(:,1) = Tzsaex(:,1)*lzsa+ lzl;
    Tlzsuex(:,2) = Tzsaex(:,2);
    Tlzmex(:,1) = Tzmex(:,1)*lzm + lza + lzsa + lzl;
    Tlzmex(:,2) = Tzmex(:,2);
else
    Tlzsuex=[];
    Tlzmex=[];
    lzsa = 0;
    lzm = 0;
end
Tzlex = read(excel,sheet1);

Tzaex = read(excel,sheet3); 

Tzvex = read(excel,sheet5); 
Tlzlex(:,1) = Tzlex(:,1)*lzl;
Tlzlex(:,2) = Tzlex(:,2);

Tlzaex(:,1) = Tzaex(:,1)*lza + lzsa + lzl;
Tlzaex(:,2) = Tzaex(:,2);

Tlzvex(:,1) = Tzvex(:,1)*lzv + lzm + lza + lzsa + lzl;
Tlzvex(:,2) = Tzvex(:,2);
data = [Tlzlex;Tlzsuex;Tlzaex;Tlzmex;Tlzvex];
[unique_positions, first_occurrence_idx] = unique(data(:, 1), 'stable');
unique_data = data(first_occurrence_idx, :);
unique_data(:,2) = unique_data(:,2);
plot(unique_data(:,1), unique_data(:,2))
function data = read(excel,sheet)
data1 = readcell(excel,'Sheet',sheet);

try
data2 = readcell(excel,'Sheet',[sheet, '<2>']);
data3 = readcell(excel,'Sheet',[sheet, '<3>']);
data4 = readcell(excel,'Sheet',[sheet, '<4>']);
catch
data2 = [];
data3 = [];
data4 = [];
end
data = [data1,data2(:,2:length(data2)),data3(:,2:length(data3)),data4(:,2:length(data4))];
data = [str2double(data(3,2:length(data))); cell2mat(data(6,2:length(data)))];
data = data';
end
end