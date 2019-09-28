f = addfile("./nodes.nc", "r")

wks = gsn_open_wks("pdf", "ngb_matrix")

dims = dimsizes(f->ngb_idx)
num_node = dims(0)
num_ngb = dims(1)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@trXMinF = 0.0
res@trXMaxF = num_node + 1.0
res@trYMinF = 0.0
res@trYMaxF = num_node + 1.0
res@trYReverse = True

res_node = True
res_node@gsMarkerColor = "blue"
res_node@gsMarkerSizeF = 1.0

base = gsn_blank_plot(wks, res)

x = new((/num_node * num_ngb/), "float")
y = new((/num_node * num_ngb/), "float")

k = 0
do i = 0, num_node - 1
  x(k:k+num_ngb-1) = f->ngb_idx(i,:)
  y(k:k+num_ngb-1) = i + 1
  k = k + num_ngb
end do

plot = gsn_add_polymarker(wks, base, x, y, res_node)

draw(base)
frame(wks)
