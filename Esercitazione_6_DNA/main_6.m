load Sample_Strings.mat

num_basi = zeros(5,4);                    %% campioni x basi
e = zeros(5,1);                           %% errori

[num_basi(1,:), e(1)] = conta_basi(S1);
[num_basi(2,:), e(2)] = conta_basi(S2);
[num_basi(3,:), e(3)] = conta_basi(S3);  %% RNA
[num_basi(4,:), e(4)] = conta_basi(S4);  %% lunghezza eccessiva
[num_basi(5,:), e(5)] = conta_basi(S5);

ris = num_basi - numB;
disp(ris)