

RawData = importdata('molecule_004000_ERI_edit.dat');


RawDataLength = size(RawData,1);
Ncont = RawData(end,1);
gabcd_ref = zeros(Ncont,Ncont,Ncont,Ncont);
for t = 1:RawDataLength
    gabcd_ref(RawData(t,1)+1,RawData(t,2)+1,RawData(t,3)+1,RawData(t,4)+1) = RawData(t,5);

end
