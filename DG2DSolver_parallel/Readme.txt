Bản backup 27/08/2020:
- Không có flux control, chỉ dùng Lax Friedrichs flux
- Bản chạy ok bậc 0 và bậc 2
- Có apply điều kiện slip

Bản backup 31/08/2020:
- Backup bản có flux control, bản này có thể sử dụng LxF, Roe và HLLE
- Bản chạy OK với tất cả các giá trị bậc chính xác
- Cập nhật hàm đọc file, các file điều kiện biên bây giờ không cần sắp xếp theo thứ tự
- Bản chạy OK với tất cả các số lượng điểm Gauss
- Có apply điều kiện slip
- Có hiệu chỉnh tại outlet để hạn chế hiện tượng backflow, tuy nhiên cần cải tiến thêm vì chưa thực sự hiệu quả
- Code được update lên github vào ngày 30/08/2020

Bản backup 1/09/2020:
- Backup bản có flux control, bản này có thể sử dụng LxF, Roe và HLLE
- Bản chạy OK với tất cả các giá trị bậc chính xác
- Cập nhật hàm đọc file, các file điều kiện biên bây giờ không cần sắp xếp theo thứ tự
- Bản chạy OK với tất cả các số lượng điểm Gauss
- Có apply điều kiện slip
- Có sửa option viscous off --> inviscid. Tăng tốc độ tính khi viscous = off
- Có hiệu chỉnh tại outlet để hạn chế hiện tượng backflow, tuy nhiên cần cải tiến thêm vì chưa thực sự hiệu quả
- Code được update lên github vào ngày 31/08/2020

Bản backup 17/09/2020:
- Sửa lại 1 số điều kiện if-else bị dư trong hàm TVDRK_1Step
- Backup để chuẩn bị dev phần mass diffusion
- Code không update lên github

Bản backup 17/11/2020:
- Sửa lại đk biên theo hướng decompose biến U -> correct biến primary theo đk biên -> reconstruct biến U- từ biến primary.
- Apply điều kiện biên zeroGradient, trong đó:
	+ Tính tang và norm grad từ phía +.
	+ ở phía trừ, tang grad - = tang grad +, norm grad - = 0
- Bị lỗi điều kiện biên -> tạm dừng.
- Code không update lên github

Bản backup 17/11/2020-lần 2:
- Update từ bản backup 17/9/2020. Sử dụng hàm zeroGradient phát triển ở bản 17/11/2020-lần 1 để correct cho div(U) tại biên outlet => xử lý được hiện tượng backflow.
- Code được update lên github vào ngày 17/11/2020