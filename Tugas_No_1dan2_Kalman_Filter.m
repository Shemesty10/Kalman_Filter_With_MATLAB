% Tugas Identifikasi dan Penyaringan
% Kalman Filter Nomor 1-2

clc
clear all

Index_k = 0:8; %Nanti harus di +1
Index_MATLAB =[0 (Index_k + 1)]; %Array Matlab harus dari 1
Posisi = [0 10 19 30 38 49 61 69 80];

Transformasi_State_A = [1 1;0 1];
Hubungan_Input_Output_H = [1 0];

Matriks_Cov_Disturb_Sistem_Q = [0 0; 0 1];
Matriks_Cov_Disturb_Pengukuran_R = 0.5;

Cov_Error_Prediksi_P_minus(:,:,1) = eye(2);

State_Prediksi_x_hat_Satu_Langkah(:,1) = [0;
                             0];

%Dapat disebut disebagai Prediksi x saat k ke 0 (array matlab 1)

for i = 1 : length(Index_k)

    Gain_Kalman_K(:,i) = Cov_Error_Prediksi_P_minus(:,:,i) * Hubungan_Input_Output_H' * inv((Hubungan_Input_Output_H * Cov_Error_Prediksi_P_minus(:,:,i) * Hubungan_Input_Output_H') + Matriks_Cov_Disturb_Pengukuran_R);
    Cov_Error_Estimasi_P(:,:,i) = (eye(2) - (Gain_Kalman_K(:,i) * Hubungan_Input_Output_H)) * Cov_Error_Prediksi_P_minus(:,:,i);
    State_Estimasi_x_hat(:,i) = State_Prediksi_x_hat_Satu_Langkah(:,i) + Gain_Kalman_K(:,i)*(Posisi(i) - (Hubungan_Input_Output_H*State_Prediksi_x_hat_Satu_Langkah(:,i)));

    n = i + 1;
    State_Prediksi_x_hat_Satu_Langkah(:,n) = Transformasi_State_A*State_Estimasi_x_hat(:,i);
    State_Prediksi_x_hat_Dua_Langkah(:,n) = (Transformasi_State_A^2)*State_Estimasi_x_hat(:,i);
    Cov_Error_Prediksi_P_minus(:,:,n) = Transformasi_State_A*Cov_Error_Estimasi_P(:,:,i)*Transformasi_State_A' + Matriks_Cov_Disturb_Sistem_Q;
    
end

State_Prediksi_x_hat_Satu_Langkah = State_Prediksi_x_hat_Satu_Langkah(:,(2:10));
State_Prediksi_x_hat_Dua_Langkah = State_Prediksi_x_hat_Dua_Langkah(:,(2:10));


%Fungsi Plot Matriks Gain Kalman K
figure(1)
Gain_Kalman_State_1 = Gain_Kalman_K(1,:);
Gain_Kalman_State_2 = Gain_Kalman_K(2,:);
plot(Index_k,Gain_Kalman_State_1)
hold on
plot(Index_k,Gain_Kalman_State_2)
grid on
xlim([min(Index_k) max(Index_k)])

xlabel('Index k')
ylabel('Gain Kalman')
legend('Gain Kalman State 1','Gain Kalman State 2')
title('Plot Gain Kalman (K)')

%Fungsi Plot Matriks Covariance Error Estimasi P
figure(2)
Variance_Error_State_1_Estimasi = squeeze(Cov_Error_Estimasi_P(1,1,:));
Variance_Error_State_2_Estimasi = squeeze(Cov_Error_Estimasi_P(2,2,:));
Covariance_Error_State_1_2_Estimasi = squeeze(Cov_Error_Estimasi_P(1,2,:));
Covariance_Error_State_2_1_Estimasi = squeeze(Cov_Error_Estimasi_P(2,1,:));

plot(Index_k,Variance_Error_State_1_Estimasi)
hold on
plot(Index_k,Variance_Error_State_2_Estimasi)
plot(Index_k,Covariance_Error_State_1_2_Estimasi)
plot(Index_k,Covariance_Error_State_2_1_Estimasi)
grid on
xlim([min(Index_k) max(Index_k)])

xlabel('Index k')
ylabel('Covariance Error Estimasi P')
legend('Variance Error State 1','Variance Error State 2','Covariance Error State 1 terhadap State 2','Covariance Error State 2 terhadap State 1')
title('Plot Covariance Error Estimasi P')

%Fungsi Plot Matriks Covariance Error Prediksi P minus
figure(3)
Variance_Error_State_1_Prediksi = squeeze(Cov_Error_Prediksi_P_minus(1,1,:));
Variance_Error_State_2_Prediksi = squeeze(Cov_Error_Prediksi_P_minus(2,2,:));
Covariance_Error_State_1_2_Prediksi = squeeze(Cov_Error_Prediksi_P_minus(1,2,:));
Covariance_Error_State_2_1_Prediksi = squeeze(Cov_Error_Prediksi_P_minus(2,1,:));

plot(Index_MATLAB,Variance_Error_State_1_Prediksi)
hold on
plot(Index_MATLAB,Variance_Error_State_2_Prediksi)
plot(Index_MATLAB,Covariance_Error_State_1_2_Prediksi);
plot(Index_MATLAB,Covariance_Error_State_2_1_Prediksi);
grid on
xlim([min(Index_MATLAB) max(Index_MATLAB)])

xlabel('Index k+1')
ylabel('Prediksi Covariance Error State P^-')
legend('Prediksi Variance Error State 1','Prediksi Variance Error State 2','Prediksi Covariance Error State 1 terhadap State 2','Prediksi Covariance Error State 2 terhadap state 1')
title('Plot Prediksi Covariance Error Estimasi P^-')

%Fungsi Plot State Posisi (Real, Estimasi, 
% Prediksi Satu Langkah, Prediksi 2 Langkah)
figure(4)

Estimasi_State_Posisi = State_Estimasi_x_hat(1,:);
Prediksi_Satu_Langkah_State_Posisi = State_Prediksi_x_hat_Satu_Langkah(1,:);
Prediksi_Dua_Langkah_State_Posisi = State_Prediksi_x_hat_Dua_Langkah(1,:);

subplot(1,2,1)
plot(Index_k,Posisi)
hold on
plot(Index_k,Estimasi_State_Posisi);
plot(Index_k+1,Prediksi_Satu_Langkah_State_Posisi);
plot(Index_k+2,Prediksi_Dua_Langkah_State_Posisi);
grid on
xlim([min(Index_k) max(Index_k)+3])

xlabel('Index k')
ylabel('Nilai State Posisi')
legend('Measured Posisi','Estimasi Posisi','Prediksi 1Langkah Posisi','Prediksi 2Langkah Posisi')
title('Pengukuran, Estimasi, Prediksi Satu Langkah, dan Prediksi Dua Langkah State Posisi')

%Fungsi Plot State Kecepatan (Estimasi, Prediksi Satu Langkah
% Prediksi 2 Langkah)
subplot(1,2,2)
Estimasi_State_Kecepatan = State_Estimasi_x_hat(2,:);
Prediksi_Satu_Langkah_State_Kecepatan = State_Prediksi_x_hat_Satu_Langkah(2,:);
Prediksi_Dua_Langkah_State_Kecepatan = State_Prediksi_x_hat_Dua_Langkah(2,:);

plot(Index_k,Estimasi_State_Kecepatan)
hold on
plot(Index_k+1,Prediksi_Satu_Langkah_State_Kecepatan);
plot(Index_k+2,Prediksi_Dua_Langkah_State_Kecepatan);
grid on
xlim([min(Index_k) max(Index_k)+3])

xlabel('Index k')
ylabel('Nilai State Kecepatan')
legend('Estimasi Kecepatan','Prediksi 1Langkah Kecepatan','Prediksi 2Langkah Kecepatan')
title('Estimasi, Prediksi Satu Langkah, dan Prediksi Dua Langkah State Kecepatan')