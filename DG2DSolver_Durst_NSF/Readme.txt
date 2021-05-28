Bản backup 22/05/2021:
- Đổi hoàn toàn sang mô hình Extended NSF Durst.
- Cần modify lại điều kiện biên vận tốc tại wall: vận tốc convection vẫn như classical NSF, vận tốc diffusion pháp tuyến với wall bằng 0.
- MÔ HÌNH DURST CHẠY OK KHÔNG BỊ CRASHHHHHHHH! REESE-GREENSHIELD LÀ HÀNG ĐỂUUUUUU.

Bản backup 27/05/2021:
- Update mô hình Durst:
	+ Dễ cài đặt vào/gỡ ra khỏi Classical NSF-DG, code của Classical NSF chỉ bị can thiệp vào ở 3 hàm.
	+ Bổ sung tích phân bị thiếu ở RHS khi giải ở bậc cao (ngoài additional term ở viscous component còn tích phân int(div(Ec)*phi) trên volume Omega. Vì div(Ec) ~ sum(Ec_i*div(phi_i)) nên khi order i=1 thì term này bằng 0) => hơi tốn não khi tìm cách tính tích phân này.
	+ Điều chỉnh điều kiện biên của vận tốc tại wall (bằng hàm bcForExtNSF_Durst::dropNormDiffVel cũng thuộc mô hình Durst) để loại bỏ thành phần vận tốc do self diffusion vuông góc với wall (đương nhiên là hiện tượng self diffusion không xảy ra trên phương vuông góc wall).
	
Bản backup 28/05/2021:
- Cải thiện hiệu suất của mô hình Durst: nhân hệ số theta2 của positivity preserving limiter vào command tính Ec_order trong hàm correctEnergyEqnVolIntTerm() vì ở bậc cao, Ec_order có thể có giá trị unphysical ở vùng có strong discontinuity -> làm solver mất ổn định.
- Sau khi cải thiện, solver có thể chạy mô hình Durst Kn=0.01 ở số Co = 0.1 (lúc trước chỉ chạy được ở Co < 0.03).