
%% Load Data
V=[];
A0=dlmread("expm-full-AdipoD0_2.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD0_2.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

Ap5=dlmread("expm-full-AdipoD0.5.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD0.5.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

A1=dlmread("expm-full-AdipoD1.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD1.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

A2=dlmread("expm-full-AdipoD2.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD2.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

A3=dlmread("expm-full-AdipoD3_2.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD3_2.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

A7=dlmread("expm-full-AdipoD7_2.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD7_2.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

A10=dlmread("expm-full-AdipoD10.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD10.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

A14=dlmread("expm-full-AdipoD14.txt",'\t',1,0);
 fid = fopen('expm-full-AdipoD14.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O0=dlmread("expm-full-OsteoD0.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD0.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

Op5=dlmread("expm-full-OsteoD0.5.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD0.5.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O1=dlmread("expm-full-OsteoD1.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD1.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O2=dlmread("expm-full-OsteoD2.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD2.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O3=dlmread("expm-full-OsteoD3.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD3.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O5=dlmread("expm-full-OsteoD5.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD5.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O7=dlmread("expm-full-OsteoD7_2.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD7_2.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O10=dlmread("expm-full-OsteoD10.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD10.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O14=dlmread("expm-full-OsteoD14.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD14.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


O21=dlmread("expm-full-OsteoD21.txt",'\t',1,0);
 fid = fopen('expm-full-OsteoD21.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];



Ch0=dlmread("expm-full-ChondroD0.txt",'\t',1,0);
 fid = fopen('expm-full-ChondroD0.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];


Ch7=dlmread("expm-full-ChondroD7.txt",'\t',1,0);
 fid = fopen('expm-full-ChondroD7.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

Ch14=dlmread("expm-full-ChondroD14.txt",'\t',1,0);
 fid = fopen('expm-full-ChondroD14.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

Ch28=dlmread("expm-full-ChondroD28.txt",'\t',1,0);
 fid = fopen('expm-full-ChondroD28.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];

Ch42=dlmread("expm-full-ChondroD42.txt",'\t',1,0);
 fid = fopen('expm-full-ChondroD42.txt');
 varNames = strsplit(fgetl(fid), '\t');
 fclose(fid);
 V=[V,varNames];
%writecell(V','MTL_Barcodes.csv')

Gene_Labels=readtable('gene-conversion.txt');
Gene_Labels=table2cell(Gene_Labels);

%% Combine data

D=[A0,Ap5,A1,A2,A3,A7,A10,A14,O0,Op5,O1,O2,O3,O5,O7,O10,O14,O21,Ch0,Ch7,Ch14,Ch28,Ch42];

DL=D;

AA0=1:length(A0(1,:));
AAp5=max(AA0)+(1:length(Ap5(1,:)));
AA1=max(AAp5)+(1:length(A1(1,:)));
AA2=max(AA1)+(1:length(A2(1,:)));
AA3=max(AA2)+(1:length(A3(1,:)));
AA7=max(AA3)+(1:length(A7(1,:)));
AA10=max(AA7)+(1:length(A10(1,:)));
AA14=max(AA10)+(1:length(A14(1,:)));
AO0=max(AA14)+(1:length(O0(1,:)));
AOp5=max(AO0)+(1:length(Op5(1,:)));
AO1=max(AOp5)+(1:length(O1(1,:)));
AO2=max(AO1)+(1:length(O2(1,:)));
AO3=max(AO2)+(1:length(O3(1,:)));
AO5=max(AO3)+(1:length(O5(1,:)));
AO7=max(AO5)+(1:length(O7(1,:)));
AO10=max(AO7)+(1:length(O10(1,:)));
AO14=max(AO10)+(1:length(O14(1,:)));
AO21=max(AO14)+(1:length(O21(1,:)));
ACh0=max(AO21)+(1:length(Ch0(1,:)));
ACh7=max(ACh0)+(1:length(Ch7(1,:)));
ACh14=max(ACh7)+(1:length(Ch14(1,:)));
ACh28=max(ACh14)+(1:length(Ch28(1,:)));
ACh42=max(ACh28)+(1:length(Ch42(1,:)));

clear A0 Ap5 A1 A2 A3 A7 A10 A14 O0 Op5 O1 O2 O3 O5 O7 O10 O14 O21 Ch0 Ch7 Ch14 Ch28 Ch42

Cat_Lab1=zeros(length(D(1,:)),1);
Cat_Lab1(AA0)=1;
Cat_Lab1(AAp5)=2;
Cat_Lab1(AA1)=3;
Cat_Lab1(AA2)=4;
Cat_Lab1(AA3)=5;
Cat_Lab1(AA7)=7;
Cat_Lab1(AA10)=8;
Cat_Lab1(AA14)=9;
Cat_Lab1(AO0)=1;
Cat_Lab1(AOp5)=2;
Cat_Lab1(AO1)=3;
Cat_Lab1(AO2)=4;
Cat_Lab1(AO3)=5;
Cat_Lab1(AO5)=6;
Cat_Lab1(AO7)=7;
Cat_Lab1(AO10)=8;
Cat_Lab1(AO14)=9;
Cat_Lab1(AO21)=10;
Cat_Lab1(ACh0)=1;
Cat_Lab1(ACh7)=3;
Cat_Lab1(ACh14)=5;
Cat_Lab1(ACh28)=7;
Cat_Lab1(ACh42)=10;
Time_Lab=categorical(Cat_Lab1,1:10,["0","0.5","1","2","3","5","7","10","14","21"]);

Cat_Lab2=zeros(length(D(1,:)),1);
Cat_Lab2(AA0)=2;
Cat_Lab2(AAp5)=2;
Cat_Lab2(AA1)=2;
Cat_Lab2(AA2)=2;
Cat_Lab2(AA3)=2;
Cat_Lab2(AA7)=2;
Cat_Lab2(AA10)=2;
Cat_Lab2(AA14)=2;
Cat_Lab2(AO0)=1;
Cat_Lab2(AOp5)=1;
Cat_Lab2(AO1)=1;
Cat_Lab2(AO2)=1;
Cat_Lab2(AO3)=1;
Cat_Lab2(AO5)=1;
Cat_Lab2(AO7)=1;
Cat_Lab2(AO10)=1;
Cat_Lab2(AO14)=1;
Cat_Lab2(AO21)=1;
Cat_Lab2(ACh0)=3;
Cat_Lab2(ACh7)=3;
Cat_Lab2(ACh14)=3;
Cat_Lab2(ACh28)=3;
Cat_Lab2(ACh42)=3;
Type_Lab=categorical(Cat_Lab2,1:3,["Osteo","Adipo","Chondro"]);

Cat_Lab3=strings(length(D(1,:)),1);
Cat_Lab3(AA0)="A0";
Cat_Lab3(AAp5)="Ap5";
Cat_Lab3(AA1)="A1";
Cat_Lab3(AA2)="A2";
Cat_Lab3(AA3)="A3";
Cat_Lab3(AA7)="A7";
Cat_Lab3(AA10)="A10";
Cat_Lab3(AA14)="A14";
Cat_Lab3(AO0)="O0";
Cat_Lab3(AOp5)="Op5";
Cat_Lab3(AO1)="O1";
Cat_Lab3(AO2)="O2";
Cat_Lab3(AO3)="O3";
Cat_Lab3(AO5)="O5";
Cat_Lab3(AO7)="O7";
Cat_Lab3(AO10)="O10";
Cat_Lab3(AO14)="O14";
Cat_Lab3(AO21)="O21";
Cat_Lab3(ACh0)="C0";
Cat_Lab3(ACh7)="C7";
Cat_Lab3(ACh14)="C14";
Cat_Lab3(ACh28)="C28";
Cat_Lab3(ACh42)="C42";

%% Remove Genes with little variability
V=var(D,0,2);
Ind=find(V>=.70);
Gene_Labels=Gene_Labels(Ind);
DD=D;


DD=DD(Ind,:);
DD=DD./sum(DD,1);

% ^^ These are the filtered genes that we care about

%% optimal number of components

%{
[u,s,v]=svd(DD);
ss=diag(s);
figure
scatter(1:length(ss)-1,diff(log(ss)))
%}

[W,H]=NMF(DD,12);

%writematrix(H([10,12,7,9,3,4,2,8,6,1,5,11],:),'MTL_Arch.csv')
%writematrix(Cat_Lab3,'MTL_Cell_Label.csv')




%% Save Archetype scores
WW=W(:,[10,12,7,9,3,4,2,8,6,1,5,11]); %reordered for biological interpretability
T=table('Size',[length(Gene_Labels),13],'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'Gene','1','2','3','4','5','6','7','8','9','10','11','12'});
T.Gene=Gene_Labels;
T{:,2:13}=WW;
%writetable(T,'Archetype_Gene_Scores.csv')

%% All Correlates with each archetype in MTL

Gene_Labels2=readtable('gene-conversion.txt');
Gene_Labels2=table2cell(Gene_Labels2);
S=std(DL,0,2);
DL=DL(S>0,:);
Gene_Labels2=Gene_Labels2(S>0);

C=corr(H',DL');
C=C';
Corr_H=table('Size',[length(Gene_Labels2) 13],'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'Gene','1','2','3','4','5','6','7','8','9','10','11','12'});
Corr_H.Gene=Gene_Labels2;
Corr_H{:,2:13}=C(:,[10,12,7,9,3,4,2,8,6,1,5,11]);

%writetable(Corr_H,'Archetype_Gene_Correlations_MTL.csv')

clear Corr_H
clear DL

%% Overlay with PDX

DL=[];
Gene_Labels_PDX=readtable('gene-conversion_pdx.txt');
Gene_Labels2=table2cell(Gene_Labels_PDX);
[~,ga,gb]=intersect(Gene_Labels,Gene_Labels2);
WW=W(ga,:);

clear AA0 AA1 AA10 AA14 AA2 AA3 AA7 AAp5 ACh0
clear ACh14 ACh28 ACh42 ACh7 AO0 AO10 AO14
clear AO1 AO2 AO21 AO3 AO5 AO7 AOp5 D DD errm

PDX1=dlmread("expm-full-PDX-1-Warm.txt",'\t',1,0);
DL=[DL,PDX1];
PDX1=PDX1(gb,:);
PDX1=PDX1./sum(PDX1,1);
PDX_H1=NMF_New_Weights(PDX1,WW,12);
PDX_H1=PDX_H1([10,12,7,9,3,4,2,8,6,1,5,11],:);

PDX2=dlmread("expm-full-PDX-2-Warm.txt",'\t',1,0);
DL=[DL,PDX2];
PDX2=PDX2(gb,:);
PDX2=PDX2./sum(PDX2,1);
PDX_H2=NMF_New_Weights(PDX2,WW,12);
PDX_H2=PDX_H2([10,12,7,9,3,4,2,8,6,1,5,11],:);

PDX3=dlmread("expm-full-PDX-3-Warm.txt",'\t',1,0);
DL=[DL,PDX3];
PDX3=PDX3(gb,:);
PDX3=PDX3./sum(PDX3,1);
PDX_H3=NMF_New_Weights(PDX3,WW,12);
PDX_H3=PDX_H3([10,12,7,9,3,4,2,8,6,1,5,11],:);


T=[PDX_H1,PDX_H2,PDX_H3];
TInd=[1*ones(1,length(PDX1(1,:))),2*ones(1,length(PDX2(1,:))),3*ones(1,length(PDX3(1,:)))];     

%writematrix(T,'PDX_Arch.csv')
%writematrix(TInd,'PDX_Cell_Label.csv')

S=std(DL,0,2);
DL=DL(S>0,:);
Gene_Labels2=Gene_Labels2(S>0);

C=corr(T',DL');
C=C';
Corr_H=table('Size',[length(Gene_Labels2) 13],'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'Gene','1','2','3','4','5','6','7','8','9','10','11','12'});
Corr_H.Gene=Gene_Labels2;
Corr_H{:,2:13}=C;

%writetable(Corr_H,'Archetype_Gene_Correlations_PDX.csv')

clear Corr_H DL


%% OS11
DL=[];
Gene_Labels_OS=readtable('gene-conversion_OS.txt');
Gene_Labels2=table2cell(Gene_Labels_OS);
[~,ga,gb]=intersect(Gene_Labels,Gene_Labels2);
WW=W(ga,:);


OS1=dlmread("expm-full-BC2.txt",'\t',1,0);
DL=[DL,OS1];
OS1=OS1(gb,:);
OS1=OS1./sum(OS1,1);
OS_H1=NMF_New_Weights(OS1,WW,12);
OS_H1=OS_H1([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS2=dlmread("expm-full-BC3.txt",'\t',1,0);
DL=[DL,OS2];
OS2=OS2(gb,:);
OS2=OS2./sum(OS2,1);
OS_H2=NMF_New_Weights(OS2,WW,12);
OS_H2=OS_H2([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS3=dlmread("expm-full-BC5.txt",'\t',1,0);
DL=[DL,OS3];
OS3=OS3(gb,:);
OS3=OS3./sum(OS3,1);
OS_H3=NMF_New_Weights(OS3,WW,12);
OS_H3=OS_H3([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS4=dlmread("expm-full-BC6.txt",'\t',1,0);
DL=[DL,OS4];
OS4=OS4(gb,:);
OS4=OS4./sum(OS4,1);
OS_H4=NMF_New_Weights(OS4,WW,12);
OS_H4=OS_H4([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS5=dlmread("expm-full-BC10.txt",'\t',1,0);
DL=[DL,OS5];
OS5=OS5(gb,:);
OS5=OS5./sum(OS5,1);
OS_H5=NMF_New_Weights(OS5,WW,12);
OS_H5=OS_H5([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS6=dlmread("expm-full-BC11.txt",'\t',1,0);
DL=[DL,OS6];
OS6=OS6(gb,:);
OS6=OS6./sum(OS6,1);
OS_H6=NMF_New_Weights(OS6,WW,12);
OS_H6=OS_H6([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS7=dlmread("expm-full-BC16.txt",'\t',1,0);
DL=[DL,OS7];
OS7=OS7(gb,:);
OS7=OS7./sum(OS7,1);
OS_H7=NMF_New_Weights(OS7,WW,12);
OS_H7=OS_H7([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS8=dlmread("expm-full-BC17.txt",'\t',1,0);
DL=[DL,OS8];
OS8=OS8(gb,:);
OS8=OS8./sum(OS8,1);
OS_H8=NMF_New_Weights(OS8,WW,12);
OS_H8=OS_H8([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS9=dlmread("expm-full-BC20.txt",'\t',1,0);
DL=[DL,OS9];
OS9=OS9(gb,:);
OS9=OS9./sum(OS9,1);
OS_H9=NMF_New_Weights(OS9,WW,12);
OS_H9=OS_H9([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS10=dlmread("expm-full-BC21.txt",'\t',1,0);
DL=[DL,OS10];
OS10=OS10(gb,:);
OS10=OS10./sum(OS10,1);
OS_H10=NMF_New_Weights(OS10,WW,12);
OS_H10=OS_H10([10,12,7,9,3,4,2,8,6,1,5,11],:);

OS11=dlmread("expm-full-BC22.txt",'\t',1,0);
DL=[DL,OS11];
OS11=OS11(gb,:);
OS11=OS11./sum(OS11,1);
OS_H11=NMF_New_Weights(OS11,WW,12);
OS_H11=OS_H11([10,12,7,9,3,4,2,8,6,1,5,11],:);


T=[OS_H1,OS_H2,OS_H3,OS_H4,OS_H5,OS_H6,OS_H7,OS_H8,OS_H9,OS_H10,OS_H11];
TInd=[1*ones(1,length(OS1(1,:))),2*ones(1,length(OS2(1,:))),3*ones(1,length(OS3(1,:))),4*ones(1,length(OS4(1,:))),5*ones(1,length(OS5(1,:))),6*ones(1,length(OS6(1,:))),7*ones(1,length(OS7(1,:))),8*ones(1,length(OS8(1,:))),9*ones(1,length(OS9(1,:))),10*ones(1,length(OS10(1,:))),11*ones(1,length(OS11(1,:)))];     

%writematrix(T,'OS_Arch.csv')
%writematrix(TInd,'OS_Cell_Label.csv')

S=std(DL,0,2);
DL=DL(S>0,:);
Gene_Labels2=Gene_Labels2(S>0);

C=corr(T',DL');
C=C';
Corr_H=table('Size',[length(Gene_Labels2) 13],'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'Gene','1','2','3','4','5','6','7','8','9','10','11','12'});
Corr_H.Gene=Gene_Labels2;
Corr_H{:,2:13}=C;

%writetable(Corr_H,'Archetype_Gene_Correlations_OS.csv')

clear Corr_H DL

%% TARGET_OS

Gene_Labels_T=readtable('gene-conversion_tos.txt');
Gene_Labels2=table2cell(Gene_Labels_T);
[~,ga,gb]=intersect(Gene_Labels,Gene_Labels2);
WW=W(ga,:);

TG=dlmread("target_os.txt",'\t',1,0);
TG=TG(gb,:);
TG=TG./sum(TG,1);
H_TG=NMF_New_Weights(TG,WW,12);
H_TG=H_TG([10,12,7,9,3,4,2,8,6,1,5,11],:);
writematrix(H_TG,'TG_Arch.csv')
