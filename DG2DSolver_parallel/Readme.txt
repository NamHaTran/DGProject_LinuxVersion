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

Bản backup 22/04/2021 (BIG UPDATE):
- Hoàn thành modify cách xử lý điều kiện biên theo hướng decompose biến U -> correct biến primary theo đk biên -> reconstruct biến U- từ biến primary.
- Tách code xử lý điều kiện biên thành folder riêng.
- Tách các hàm thuộc về phần parallel processing ra thành 1 folder.
- Chuyển phương thức send/recv data từ send/recv theo cell data sang send/recv theo điểm Gauss trên các edge biên match => giảm thời gian tính toán vì không cần phải tính lại giá trị điểm Gauss thuộc cell neighbor (thuộc processor neighbor) ở processor đang xét. Đồng thời giải quyết bug missmatch khi giải T ở mode implicit (mass diffusion = true).
- Tách vector base Gauss points và base Gauss weight thành 2 cho mỗi vector (1 cho volume Gauss, 1 cho surface Gauss), tuy nhiên chưa tách hoàn toàn volume và surface Gauss ra để dùng theo 2 setting số điểm Gauss riêng biệt. Bước hiện tại là tiền đề để làm việc này => giúp tăng tốc độ tính (giảm số điểm Gauss trong khi vẫn đảm bảo bậc chính xác)
- Update ghi chú Doxygen.

Bản backup 09/05/2021:
- Update hàm mass difusion limiter: về vơ bản đã chạy được, hiện tại đang trong giai đoạn test và tối ưu hóa code. Cần phải test thêm về tác dụng của limiter này trong trường hợp Kn - 0.05 và 0.1 (đã chạy tốt với Kn = 0.01, có thể chạy với max Co = 0.01, lúc trước 0.001 vẫn gây crash).
- Bổ sung một số update nhỏ ở phần hiển thị limiter setting.
- Code được update lên github vào ngày 09/05/2021

Bản backup 10/05/2021:
- Update hàm đọc limiter setting để đọc mass diffuson setting.
- Update hàm show setting của limiter.
- Modify hàm tính pointAuxValue để apply limiter vào cả 4 biến Sm = div(U) khi mass diffusion = ON.
- Update code lên github cùng ngày.