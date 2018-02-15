%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Dichromatic reflectance model %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read the observations with large pixels 
A = imread('kobi.png');
figure(1);
imshow(A);
%% Use the superpixel to segment the picture
[L,N] = superpixels(A,1000);
observations = zeros(N,3);
outputImage = zeros(size(A),'like',A);
idx = label2idx(L);
numRows = size(A,1);
numCols = size(A,2);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = mean(A(redIdx));
    outputImage(greenIdx) = mean(A(greenIdx));
    outputImage(blueIdx) = mean(A(blueIdx));
    observations(labelVal,:) = [mean(A(redIdx)) mean(A(greenIdx)) mean(A(blueIdx))];
end    
figure(2);
imshow(outputImage,'InitialMagnification',67);
%% Use the Gauss-Seidel iterations to minimize the cost function
% Initialize the parameters here
C_l = observations;
C_d = 100*ones(N,3);
M_d = ones(1,N);
C_s = 100*ones(1,3);
M_s = ones(1,N);
for p = 1:1000
%         for j = 1:N
%         M_d(j) = C_d(j,:)*(C_l(j,:) - M_s(j)*C_s)'/(C_d(j,:)*C_d(j,:)');
%         end
% 
%         for j = 1:N
%         M_s(j) = C_s*(C_l(j,:)-M_d(j)*C_d(j,:))'/(C_s*C_s');
%         end 
% 
%     C_s(1) = sum(M_s.*(C_l(:,1)'-M_d.*C_d(:,1)'))/(sum(M_s.^2));    
%     C_s(2) = sum(M_s.*(C_l(:,2)'-M_d.*C_d(:,2)'))/(sum(M_s.^2));   
%     C_s(3) = sum(M_s.*(C_l(:,3)'-M_d.*C_d(:,3)'))/(sum(M_s.^2));   
%     
%     for j = 1:N
%     C_d(j,1) = M_d(j)*(C_l(j,1) - M_s(j)*C_s(1))/M_d(j)^2; 
%     C_d(j,2) = M_d(j)*(C_l(j,2) - M_s(j)*C_s(2))/M_d(j)^2;
%     C_d(j,3) = M_d(j)*(C_l(j,3) - M_s(j)*C_s(3))/M_d(j)^2;
%     end

        for j = 1:N
        M_d(j) = C_d(j,:)*(M_s(j)*C_s-C_l(j,:))'/(C_d(j,:)*C_d(j,:)');
        end

        for j = 1:N
        M_s(j) = C_s*(M_d(j)*C_d(j,:)-C_l(j,:))'/(C_s*C_s');
        end 

    C_s(1) = sum(M_s.*(M_d.*C_d(:,1)'- C_l(:,1)'))/(sum(M_s.^2));    
    C_s(2) = sum(M_s.*(M_d.*C_d(:,2)'- C_l(:,2)'))/(sum(M_s.^2));   
    C_s(3) = sum(M_s.*(M_d.*C_d(:,3)'- C_l(:,3)'))/(sum(M_s.^2));   
    
    for j = 1:N
    C_d(j,1) = M_d(j)*(M_s(j)*C_s(1) - C_l(j,1))/M_d(j)^2; 
    C_d(j,2) = M_d(j)*(M_s(j)*C_s(2) - C_l(j,2))/M_d(j)^2;
    C_d(j,3) = M_d(j)*(M_s(j)*C_s(3) - C_l(j,3))/M_d(j)^2;
    end

end


%% Output the true color of objects
TrueImage = zeros(size(A),'like',A);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    TrueImage(redIdx) = C_d(labelVal,1);
    TrueImage(greenIdx) = C_d(labelVal,2);
    TrueImage(blueIdx) = C_d(labelVal,3);
end   
figure(3);
imshow(TrueImage,'InitialMagnification',67);

