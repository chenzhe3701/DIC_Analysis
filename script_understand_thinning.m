A=A;
B=(bwmorph(A,'thin',1));
C=(bwmorph(A,'thin',2));
D=(bwmorph(A,'thin',3));
E=(bwmorph(A,'thin',4));
F=(bwmorph(A,'thin',5));
figure;imagesc(A+B+C+D+E+F);

