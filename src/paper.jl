function paper()
  pyplot()
  @info "==== PLOTTING SOLITONS ===="
  plot_solitons()
  @info "==== PLOTTING TILES    ===="
  plot_tiles()
  @info "==== PLOTTING LINES    ===="
  view_all_lines()
  @info "==== PLOTTING CHEMPOT  ===="
  compare_chempot()
  nothing
end