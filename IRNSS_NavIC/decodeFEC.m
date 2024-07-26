function navData = decodeFEC(data)

% Ghi theo cot doc theo hang 73*8
for i = 1:length(data)/600
    % Tách từ phần tử 17 đến 600
    subset = data(600*(i-1)+17 : 600*i);

    % Chuyển subset thành ma trận 73x8 theo cột
    reshapedMatrix = reshape(subset, 73, 8);

    % Đọc lại ma trận theo hàng
    readByRow = reshape(reshapedMatrix', 1, []);

    data(600*(i-1)+17 : 600*i) = readByRow(1:584);
end

% Tim du lieu goc
navData = [];
for i = 1:length(data)/600
    bits = data(600*(i-1)+17 : 600*i);

    % Xác định cấu trúc trellis từ các thông số mã hóa
    constraintLength = 7; % Độ dài ràng buộc
    poly1 = [171]; % Đa thức sinh G1
    poly2 = [133]; % Đa thức sinh G2
    trellis = poly2trellis(constraintLength, [poly1, poly2]);

    % Dãy bit đã được mã hóa (ví dụ)
    receivedData = bits;

    % Giải mã đoạn bit sử dụng thuật toán Viterbi
    decodedData = vitdec(receivedData, trellis, 7, 'term', 'hard');

    if all(decodedData(end-5:end) == 0)
        navData = horzcat(navData, decodedData);
    end

end


