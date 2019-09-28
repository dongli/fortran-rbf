f = addfile("./nodes.nc", "r")

query_node = 1000

wks = gsn_open_wks("pdf", "nodes")

res = True
res@gsnDraw = False
res@gsnFrame = False
res@trXMinF = 0.0
res@trXMaxF = 1.0
res@trYMinF = 0.0
res@trYMaxF = 1.0

res_node = True
res_node@gsMarkerColor = "blue"
res_node@gsMarkerSizeF = 1.0

res_ngb = True
res_ngb@gsMarkerColor = "red"
res_ngb@gsMarkerSizeF = 1.0
res_ngb@gsMarkerIndex = 4

res_query = True
res_query@gsMarkerSizeF = 1.0
res_query@gsMarkerIndex = 12

base = gsn_blank_plot(wks, res)

plot_node = gsn_add_polymarker(wks, base, f->x(:,0), f->x(:,1), res_node)

plot_ngb = gsn_add_polymarker(wks, base, f->x(f->ngb_idx(query_node,:)-1,0), f->x(f->ngb_idx(query_node,:)-1,1), res_ngb)

plot_query = gsn_add_polymarker(wks, base, f->x(query_node,0), f->x(query_node,1), res_query)

draw(base)
frame(wks)
