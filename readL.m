function [lzl,lzsa,lza,lzm,lzv] = readL(excel,sheet1,sheet3,sheet5,varargin)
if length(varargin) == 2
sheet2 = varargin{1};
sheet4 = varargin{2};
lzsadata = readtable(excel,'Sheet',sheet2);
lzmdata = readtable(excel,'Sheet',sheet4);
lzsa = table2array(lzsadata(2,2));
lzm = table2array(lzmdata(2,2));
else
    lzsa = 0;
    lzm = 0;
end
lzldata = readtable(excel,'Sheet',sheet1);

lzadata = readtable(excel,'Sheet',sheet3);

lzvdata = readtable(excel,'Sheet',sheet5);
lzl = table2array(lzldata(2,2));

lza = table2array(lzadata(2,2));

lzv = table2array(lzvdata(2,2));
end