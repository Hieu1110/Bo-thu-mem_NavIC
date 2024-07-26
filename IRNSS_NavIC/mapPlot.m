function mapPlot(navSolutions)

% Tạo một bản đồ webmap
webmap('World Imagery');

% Thêm điểm trên bản đồ
lat1 = mean(navSolutions.latitude); % latitude trung bình của các vị trí tính toán
lon1 = mean(navSolutions.longitude); % longitude
wmPoint1 = wmmarker(lat1, lon1, 'Color', 'blue', 'FeatureName', 'Vị trí tính toán', 'OverlayName', 'Vị trí tính toán', 'IconScale', 1.5);

% Thêm điểm thứ hai trên bản đồ
lat2 = 21.00452250; % Latitude của vị trí 2
lon2 = 105.84634683; % Longitude của vị trí 2
wmPoint2 = wmmarker(lat2, lon2, 'Color', 'red', 'FeatureName', 'Vị trí thực tế', 'OverlayName', 'Vị trí thực tế', 'IconScale', 1.5);

% Sử dụng hàm 'distance' của MATLAB để tính khoảng cách
distanceInDegrees = distance(lat1, lon1, lat2, lon2); 
% Khoảng cách trả về là độ, cần chuyển đổi thành km (1 độ tương đương với khoảng 111 km)
distanceInKm = deg2km(distanceInDegrees);   

% Hiển thị khoảng cách
disp(['Khoảng cách giữa hai vị trí là: ', num2str(distanceInKm*1000), ' m']);

% Nối hai điểm với nhau trên bản đồ
wmLine = wmline([lat1, lat2], [lon1, lon2], 'FeatureName', num2str(distanceInKm*1000) + " m", 'OverlayName', "Khoảng cách: " + num2str(distanceInKm*1000) + " m", 'LineWidth', 5);

end

